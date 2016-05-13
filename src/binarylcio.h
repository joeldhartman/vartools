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
#ifndef __BINLCIO_H
#define __BINLCIO_H

enum outtype {INT=0, PACKED_FLOAT=1, UNPACKED_FLOAT=2};

typedef struct {
	enum outtype tp;
	void *nan_repr;
	unsigned char has_sign;
	unsigned short numbits;
        int bytesfromstart;
	long offset;
	double precision;
}  BinLC_OutputFormat;

/*extern struct OutputFormat *lc_column_format;
extern unsigned lc_num_columns; 
extern unsigned lc_record_size;*/

/* Writes an unsigned integer to the first num_bits bits of dest. Currently
 * assumes num_bits is divisible by 8. */
void push_unsigned_int(char *dest, int num_bits, unsigned long val);

/* Writes a signed integer to the first num_bits bits of dest. Currently 
 * assumes num_bits is divisible by 8. */
void push_signed_int(char *dest, int num_bits, long val);

/* Writes a floating point value to the first num_bits bits of dest.
 * Currently only works if num_bits is the size of either float or double. */
void push_float(char *dest, int num_bits, const double val);

/* Converts the given floating point value to integer by first dividing it
 * by the given precision, rounding and sibstracting the given offset, and
 * writes the result to dest. Currently only supporst num_bits values
 * divisible by 8. */
void pack_float(char *dest, int num_bits, const double precision, 
		const long offset, int has_sign, const double val, char *nan_repr);

/* Returns an unsigned integer derived from the first num_bits bits of 
 * source. Currently assumes num_bits is divisible by 8. */
unsigned long pop_unsigned_int(const char *source, int num_bits);

/* Returns a pointer to a signed integer derived from the first num_bits
 * bits of source. Currently assumes num_bits is divisible by 8. */
long pop_signed_int(char *source, int num_bits);

/* Reads a floating point value from source. 
 * Currently only works if num_bits is the size of either float or double. */
double pop_float(char *source, int num_bits);

/* Reads a packed floating point value from source. Currently only supports
 * num_bits values divisible by 8. */
double unpack_float(char *source, int num_bits, const double precision, 
		const long offset, int has_sign, void *nan_repr);

/* Returns a representation of nan suitable for a numbits packed floating
 * point value if the underlying integer representation is signed. Currently
 * only works for numbits divisible by 8.*/
void *signed_nan(unsigned short numbits);

/* Returns a representation of nan suitable for a numbits packed floating
 * point value if the underlying integer representation is not signed. 
 * Currently only works for numbits divisible by 8.*/
void *unsigned_nan(unsigned short numbits);

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
		char *error_message);

/* Sets one character past the end of the keyword to '\0', valstart to the
 * first character of the value, one past the last character of the
 * value to '\0' and similarly for comm_(start/end) of the given
 * header record formatted as value = keyword / comment with optional space
 * padding around the = and / characters and after the end of the comment. If
 * the given header record is empty all the keywords are set to NULL. 
 * Also, for non-empty records sets the size of the keyword (the number of
 * characters before =). */
void bin_lc_parse_header_record(char *hdr_record, int record_len, 
		char **value, char **comment, int *kwsize);

/* Releases all memory held by the current format and allocates space for a
 * new one. */
/*void reset_lc_format(unsigned numcol);*/

/* Releases the memory held by the lc_column_format array. */
void free_lc_format(BinLC_OutputFormat *lc_column_format);

/* Calculates how many characters should be skipped over in order to jump
 * over all columns between column1 and column2 not including either column1
 * or column2. */
long bits_to_skip_between(int column1, int column2, BinLC_OutputFormat *lc_column_format);

#endif
