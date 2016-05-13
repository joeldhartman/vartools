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
#ifndef DEF_LINKEDLIST
#define DEF_LINKEDLIST

struct ListNode {
  int val;
  struct ListNode *priornode;
  struct ListNode *nextnode;
};

struct LinkedList {
  int Nnodes;
  struct ListNode *firstnode;
  struct ListNode *lastnode;
};

void InitLinkedList(struct LinkedList *list);
void PushNode(struct LinkedList *list, int val);
int PopNode(struct LinkedList *list);
void RemoveNode(struct LinkedList *list, struct ListNode *node);
#endif
