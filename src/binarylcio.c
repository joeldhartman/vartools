/*     This file is part of VARTOOLS version 1.31                      */
/*                                                                           */
/*     VARTOOLS is free software: you can redistribute it and/or modify      */
/*     it under the terms of the GNU General Public License as published by  */
/*     the Free Software Foundation, either version 3 of the License, or     */
/*     (at your option) any later version.                                   */
/*                                                                           */
/*     This program is distributed in the hope that it will be useful,       */
/*     but WITHOUT ANY WARRANTY; without even the implied warranty of        */
/*     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         */
/*     GNU General Public License for more details.                          */
/*                                                                           */
/*     You should have received a copy of the GNU General Public License     */
/*     along with this program.  If not, see <http://www.gnu.org/licenses/>. */
/*                                                                           */
/*     Copyright 2007, 2008, 2009  Joel Hartman                              */
/*                                                                           */
/* This file contains code to support a binary light curve format 
   developed by Kaloyan Penev. The code in this file was written
   by Kaloyan Penev for the HATSouth project, with
   modifications by Joel Hartman to make it compatible with VARTOOLS.
*/
#include "binarylcio.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

//struct OutputFormat *lc_column_format=NULL;
//unsigned lc_num_columns=0, lc_record_size=0;

/* Writes an unsigned integer to the first num_bits bits of dest. Currently
 * assumes num_bits is divisible by 8. */
void push_unsigned_int(char *dest, int num_bits, unsigned long val)
{
	int i;
	char *source=(char*)(&val);
	for(i=0; i<num_bits/8; i++) dest[i]=source[i];
}

/* Writes a signed integer to the first num_bits bits of dest. Currently 
 * assumes num_bits is divisible by 8. */
void push_signed_int(char *dest, int num_bits, long val)
{
	if(val<0) val=(-val|1<<(num_bits-1));
	push_unsigned_int(dest, num_bits, val);
}

/* Writes a floating point value to the first num_bits bits of dest.
 * Currently only works if num_bits is the size of either float or double. */
void push_float(char *dest, int num_bits, const double val)
{
	int i, size=num_bits/8;
	char *source;
	float float_val=val;
	if(size==sizeof(double)) source=(char*)(&val);
	else source=(char*)(&float_val);
	for(i=0; i<size; i++) dest[i]=source[i];
}

/* Converts the given floating point value to integer by first dividing it
 * by the given precision, rounding and sibstracting the given offset, and
 * writes the result to dest. Currently only supporst num_bits values
 * divisible by 8. */
void pack_float(char *dest, int num_bits, const double precision, 
		const long offset, int has_sign, const double val, char *nan_repr)
{
	if(nan_repr && (isnan(val) || isinf(val))) {
		if(has_sign) push_signed_int(dest, num_bits, *(long *)(nan_repr));
		else push_unsigned_int(dest, num_bits, *(unsigned long *)(nan_repr));
	} else {
		if(has_sign) push_signed_int(dest, num_bits, 
				(long)(floor(0.5+val/precision))-offset);
		else push_unsigned_int(dest, num_bits, 
				(unsigned long)(floor(0.5+val/precision))-offset);
	}
}

/* Returns an unsigned integer derived from the first num_bits bits of 
 * source. Currently assumes num_bits is divisible by 8. */
unsigned long pop_unsigned_int(const char *source, int num_bits)
{
	int long_size=sizeof(unsigned long), i;
	unsigned long result=*(unsigned long *)(source);
	char *result_bytes=(char*)(&result);
	for(i=num_bits/8; i<long_size; i++) result_bytes[i]=0;
	return result;
}

/* Returns a pointer to a signed integer derived from the first num_bits
 * bits of source. Currently assumes num_bits is divisible by 8. */
long pop_signed_int(char *source, int num_bits)
{
	unsigned long result=pop_unsigned_int(source, num_bits);
	unsigned long test_bit=(unsigned long)(1)<<(num_bits-1);
	if (result & test_bit) {
		result^=test_bit;
		if(!result) return -test_bit;
		else return -result;
	} else return result;
}

/* Reads a floating point value from source. 
 * Currently only works if num_bits is the size of either float or double. */
double pop_float(char *source, int num_bits)
{
	int i, size=num_bits/8;
	float float_val;
	double double_val;
	char *dest;

	if(size==sizeof(double)) dest=(char*)(&double_val);
	else dest=(char*)(&float_val);
	for(i=0; i<size; i++) dest[i]=source[i];
	if(size==sizeof(double)) return double_val;
	else return float_val;
}

/* Reads a packed floating point value from source. Currently only supports
 * num_bits values divisible by 8. */
double unpack_float(char *source, int num_bits, const double precision, 
		const long offset, int has_sign, void *nan_repr)
{
	long sval;
	unsigned long uval;

	if(has_sign) {
		sval=pop_signed_int(source, num_bits);
		if(nan_repr && sval==*(long*)(nan_repr)) return 0.0/0.0;
		else return (sval+offset)*precision;
	} else {
		uval=pop_unsigned_int(source, num_bits);
		if(nan_repr && uval==*(unsigned long*)(nan_repr)) return 0.0/0.0;
		else return (uval+offset)*precision;
	}
}

/* Returns a representation of nan suitable for a numbits packed floating
 * point value if the underlying integer representation is signed. Currently
 * only works for numbits divisible by 8.*/
void *signed_nan(unsigned short numbits)
{
	void *result=malloc(sizeof(long));
	*(long*)(result)=-(unsigned long)(1)<<(numbits-1);
	return result;
}

/* Returns a representation of nan suitable for a numbits packed floating
 * point value if the underlying integer representation is not signed. 
 * Currently only works for numbits divisible by 8.*/
void *unsigned_nan(unsigned short numbits)
{
	char *result=calloc(1, sizeof(unsigned long));
	int i;

	for(i=0; i<numbits/8; i++) result[i]=~'\0';
	return result;
}

/* Parses a column format string, filling up the members of the format
 * argument :
 * - FP,## - ## bit floating point
 * - PACKED,[IF],[TF],[+-],###,#(11),### - integer packed column where the 
 *   comma separated entries in order of appearance are:
 *   - a single character, either I if the column contains integer numbers or
 *     F if it containts floating point numbers
 *   - a single character, either T if the column is allowed to contain 
 *     missing/undefined values or F if it is not
 *   - a single character, either + if the packed values correspond to 
 *     unsigned integers or - if they have a sign
 *   - exactly 3 digit unsigned integer (0 left padded) giving the number of
 *     bits used to represent the (scaled), shifted, (rounded) values
 *   - exactly 11 character long integer giving the offset applied before 
 *     packing the column but after rescaling
 *   - exactly 3 character long signed integer specifying the number of 
 *     digits preserved, must begin with a '+' or a '-'. 
 * - INT,## - a shortcut for: PACKED,I,F,-,##,0          ,+0
 * Returns 0 on success and nonzero on failure, setting the appropriate
 * python exception. */
int parse_format_string(const char *fmt_string, BinLC_OutputFormat *format,
		char *error_message)
{
	char packed_typec, nanc, signc;
	int preserve_digits;

	if(strncmp(fmt_string, "FP,", 3)==0) {
		format->tp=UNPACKED_FLOAT;
		format->numbits=atoi(fmt_string+3);
		format->nan_repr=NULL;
		if(format->numbits/8!=sizeof(double) && 
				format->numbits/8!=sizeof(float)) {
			if(error_message)
				sprintf(error_message, 
						"Unpacked floating point column format '%s' has bit "
						"size (%d) that is neither that of float (%ld) nor "
						"that of double (%ld)", fmt_string, format->numbits, 
						sizeof(float)*8, sizeof(double)*8);
			return 1;
		}
	} else if(strncmp(fmt_string, "INT,", 4)==0) {
		format->tp=INT;
		format->numbits=atoi(fmt_string+4);
		format->nan_repr=NULL;
		format->has_sign=1;
		format->offset=0;
		format->precision=1;
	} else if(strncmp(fmt_string, "PACKED,", 7)==0) {
		if(sscanf(fmt_string+7, "%c,%c,%c,%03hu,%011ld,%3d", &packed_typec, 
					&nanc, &signc, &format->numbits, &format->offset, 
					&preserve_digits)!=6) {
			if(error_message)
				sprintf(error_message, "Invalid PACKED data type in column "
						"format string: '%s'. Should be of the form: "
						"PACKED,[IF],[TF],[+-],###,#(11),###.", fmt_string);
			return 1;
		}
		format->precision=pow(10, -preserve_digits);
		switch(signc) {
			case '+': format->has_sign=0; break;
			case '-': format->has_sign=1; break;
			default: if(error_message)
						 sprintf(error_message, "Invalid sign flag: '%c' in "
								 "column format string: '%s'. Should be "
								 "either '+' or '-'.", signc, fmt_string);
					 return 1;
		}
		switch(packed_typec) {
			case 'I': format->tp=INT; format->nan_repr=NULL; break;
			case 'F': format->tp=PACKED_FLOAT; 
					  switch(nanc) {
						  case 'T': 
							  format->nan_repr=(signc=='+' ? 
									  unsigned_nan(format->numbits) :
									  signed_nan(format->numbits));
							  break;
						  case 'F': format->nan_repr=NULL; break;
						  default: if(error_message) sprintf(error_message, 
										   "Invalid nan flag: '%c' in column"
										   " format string: '%s'. Should be "
										   "either 'T' or 'F'.", nanc, 
										   fmt_string); 
									   return 1;
					  }; 
					  break;
			default: if(error_message) sprintf(error_message, 
							 "Invalid packed data type: '%c' in column "
							 "format string: '%s'. Should be either 'I' or "
							 "'F'.", packed_typec, fmt_string);
					 return 1;
		}
	} else {
		if(error_message) sprintf(error_message, 
				"Invalid data type: in column format string: '%s'. Should be"
				" one of FP, INT, PACKED.", fmt_string);
		return 1;
	}
	if(format->numbits%8) {
		if(error_message) sprintf(error_message, 
				"Bit size of column format '%s' not divisible by 8: %d", 
				fmt_string, format->numbits);
		return 1;
	}
	return 0;
}


/* Sets one character past the end of the keyword to '\0', valstart to the
 * first character of the value, one past the last character of the
 * value to '\0' and similarly for comm_(start/end) of the given
 * header record formatted as value = keyword / comment with optional space
 * padding around the = and / characters and after the end of the comment. If
 * the given header record is empty all the keywords are set to NULL. 
 * Also, for non-empty records sets the size of the keyword (the number of
 * characters before =). */
void bin_lc_parse_header_record(char *hdr_record, int record_len, 
		char **value, char **comment, int *kwsize)
{
	char *pos;
	*value=strchr(hdr_record, '=');
	if(*value==NULL) {
		*comment=NULL;
		return;
	}
	*kwsize=*value-hdr_record;
	pos=*value;
	(*value)++;
	while(**value==' ') (*value)++;
	while(pos[-1]==' ') pos--;
	*pos='\0';
	*comment=strchr(*value, '/');
	if(**value=='/') **value='\0';
	else if(*comment) {
		pos=*comment;
		while(pos[-1]==' ') pos--;
		*pos='\0';
		(*comment)++;
		while(**comment==' ') (*comment)++;
		if(**comment) {
			pos=hdr_record+record_len-1;
			while(pos[-1]==' ') pos--;
			*pos='\0';
		}
	} else {
		pos=hdr_record+record_len-1;
		while(pos[-1]==' ') pos--;
		*pos='\0';
		*comment=pos;
	}
}

/* Releases the memory held by the lc_column_format array. */
void free_lc_format(BinLC_OutputFormat *lc_column_format)
{
	unsigned i;
	free(lc_column_format->nan_repr);
}

/* Releases all memory held by the current format and allocates space for a
 * new one. */
/*
void reset_lc_format(unsigned numcol, BinLC_OutputFormat *lc_column_format)
{
	int i;

	if((*lc_column_format)) free_lc_format(lc_column_format);
	(*lc_num_columns)=numcol;
	lc_record_size=0;
	lc_column_format=(struct OutputFormat *)malloc(
			numcol*sizeof(struct OutputFormat));
	for(i=0; i<numcol; i++) lc_column_format[i].nan_repr=NULL;
}
*/

/* Calculates how many characters should be skipped over in order to jump
 * over all columns between column1 and column2 not including either column1
 * or column2. */
long bits_to_skip_between(int column1, int column2, 
			  BinLC_OutputFormat *lc_column_format)
{
	long jump=0;
	for(column1++;column1<column2; column1++) 
		jump+=lc_column_format[column1].numbits;
	return jump;
}
