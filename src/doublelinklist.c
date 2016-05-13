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
#include <stdlib.h>
#include "doublelinklist.h"

void InitLinkedList(struct LinkedList *list) {
  list->Nnodes = 0;
  list->firstnode = NULL;
  list->lastnode = NULL;
}

void PushNode(struct LinkedList *list, int val) {
  struct ListNode *tmp;
  tmp = (struct ListNode *) malloc(sizeof(struct ListNode));
  tmp->val = val;
  tmp->nextnode = NULL;
  if(!list->Nnodes) {
    list->firstnode = tmp;
  }
  tmp->priornode = list->lastnode;
  if(tmp->priornode != NULL) {
    tmp->priornode->nextnode = tmp;
  }
  list->lastnode = tmp;
  list->Nnodes += 1;
}

int PopNode(struct LinkedList *list) {
  int retval;
  struct ListNode *tmp;
  if(list->Nnodes == 0) {
    return -1;
  } else {
    tmp = list->lastnode;
    retval = list->lastnode->val;
    list->lastnode = tmp->priornode;
    if(list->lastnode != NULL) {
      list->lastnode->nextnode = NULL;
    }
    free(tmp);
    list->Nnodes -= 1;
    return retval;
  }
}

void RemoveNode(struct LinkedList *list, struct ListNode *node) {
  struct ListNode *tmp1;
  struct ListNode *tmp2;
  if(node == NULL) return;
  tmp1 = node->priornode;
  tmp2 = node->nextnode;

  if(tmp1 != NULL) {
    tmp1->nextnode = tmp2;
  }
  if(tmp2 != NULL) {
    tmp2->priornode = tmp1;
  }

  if(list->lastnode == node) {
    list->lastnode = tmp1;
  }
  if(list->firstnode == node) {
    list->firstnode = tmp2;
  }
  list->Nnodes = list->Nnodes - 1;
  free(node);
}

    
