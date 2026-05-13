//#######################################################################
//###########   READ XML ODB DUMP FROM START OF DATA FILE   #############
//#######################################################################
// the odb dump has already been read into memory
// run through the data, saving anything of interest
// *there is no error recovery from corrupted data*
//  (no bad odb dumps have been seen so far)
// 
// NOTE xml files *can* contain ">" in data - These routines fail in that case
//    [usually "&gt;" is used instead of ">", but this is not mandatory]
//    [but midas always uses "&gt;", so these routines are good for midas]
//
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>
#include "odb.h"

static char *read_xml_data(char *xml_file, Xml_node *current, char *name_ptr, int namelen, char *attr, char *attr_end);
static int free_odb_attrlist(Xml_attr *current);
static Xml_node *alloc_odb_element(char *name, int namelen);
static int free_odb_element(Xml_node *current);
static int free_odb_tree(Xml_node *root);
static int insert_odb_attr(Xml_node *current, char *name_ptr, int namelen, char *value, int vlen);
static Xml_attr *get_attr(Xml_attr *attr, char *id);
static Xml_node *get_odb_node(char *path);
static char *attr_value(Xml_attr *attr);
static char *node_value(Xml_node *node);
static int arraystr_idx(char *str, char *array, int size, int num);

//extern char midas_runtitle[SYS_PATH_LENGTH];
Xml_node *odb_tree;

#define ODB_DATATYPES  16
#define DTYPE_STRLEN   16
char odb_typestr[ODB_DATATYPES][DTYPE_STRLEN] = {
   "BYTE", "SBYTE", "CHAR", "WORD",       "SHORT", "DWORD", "INT", "BOOL",
   "FLOAT", "DOUBLE", "BITS", "STRING",    "ARRAY", "STRUCT", "KEY", "LINK"
};

// read odb value into integer variable (even if source type is float)
static Xml_val *get_odb_value(char *path, int *type);
int odbval_int(char *key, int *dst)
{
   int type;  float f_val; Xml_val *val;

   if( (val = get_odb_value(key, &type)) == NULL ){ return(-1); }
   if( type >= 0 && type <= 7 ){ // ints
      if( sscanf(val->p_val, "%d", dst) >= 1 ){ return(0); }
      printf("odbval_int:cant read integer value[%s]\n", val->p_val);
      return(-1);
   }
   if( type == 8 || type == 9 ){ // floats
      if( sscanf(val->p_val, "%f", &f_val) >= 1 ){ *dst = f_val; return(0); }
      printf("odbval_int:cant read value[%s]\n", val->p_val);
      return(-1);
   }
   printf("odbval_int[%s] Non-Numeric key-type [%d]\n", key, type);
   return(-1);
}
// read odb value into floating-point variable (even if source type is integer)
int odbval_float(char *key, float *dst)
{
   int type, i_val;  Xml_val *val;

   if( (val = get_odb_value(key, &type)) == NULL ){ return(-1); }
   if( type >= 0 && type <= 7 ){ // ints
      if( sscanf(val->p_val, "%d", &i_val) >= 1 ){ *dst = i_val; return(0); }
      printf("odbval_int:cant read integer value[%s]\n", val->p_val);
      return(-1);
   }
   if( type == 8 || type == 9 ){ // floats
      if( sscanf(val->p_val, "%f", dst) >= 1 ){ return(0); }
      printf("odbval_int:cant read value[%s]\n", val->p_val);
      return(-1);
   }
   printf("odbval_float[%s] Non-Numeric key-type [%d]\n", key, type);
   return(-1);
}

int odbarray_int(char *key, int *dst, int maxlen)
{
   int i, *i_array, type, nval, size;
   void *array;  float *f_array;

   if( get_odb_array(key, &array, &type, &nval, &size) == -1 ){ return(-1); }
   if( type >= 0 && type <= 7 ){ // integer values
      i_array = array;
      for(i=0; i<nval && i<maxlen; i++){ dst[i] = i_array[i]; }
      if( nval > maxlen ){
         printf("odbarray_int[%s] Too many entries %d [max %d]\n",
                                                  key, nval, maxlen);
         return(-1);
      }
   } else if ( type == 8 || type == 9 ){ // float values
      f_array = array;
      for(i=0; i<nval && i<maxlen; i++){ dst[i] = f_array[i]; }
      if( nval > maxlen ){
         printf("odbarray_int[%s] Too many entries %d [max %d]\n",
                                                  key, nval, maxlen);
         return(-1);
      }
   } else {
      printf("odbarray_int[%s] Non-Numeric key-type [%d]\n", key, type);
      return(-1);
   }
   return(0);
}

int odbarray_float(char *key, float *dst, int maxlen)
{
   int i, *i_array, type, nval, size;
   void *array;  float *f_array;

   if( get_odb_array(key, &array, &type, &nval, &size) == -1 ){ return(-1); }
   if( type >= 0 && type <= 7 ){ // integer values
      i_array = array;
      for(i=0; i<nval && i<maxlen; i++){ dst[i] = i_array[i]; }
      if( nval > maxlen ){
         printf("odbarray_float[%s] Too many entries %d [max %d]\n",
                                                  key, nval, maxlen);
         return(-1);
      }
   } else if ( type == 8 || type == 9 ){ // float values
      f_array = array;
      for(i=0; i<nval && i<maxlen; i++){ dst[i] = f_array[i]; }
      if( nval > maxlen ){
         printf("odbarray_float[%s] Too many entries %d [max %d]\n",
                                                  key, nval, maxlen);
         return(-1);
      }
   } else {
      printf("odbodbarray_float[%s] Non-Numeric key-type [%d]\n", key, type);
      return(-1);
   }
   return(0);
}

#define MAX_ARRAYSTR_SIZE 256
// values are in subtree, each value is own key, with index specified
// alloc resulting array here - free in caller
int get_odb_array(char *path, void **value, int *type, int *nvalues, int *size)
{
   int i, i_val, index, *i_array;
   float f_val, *f_array;
   char *tmp, *result;
   Xml_node *node;
   Xml_attr *attr;
   *size = 4;

   // find array node, and read attrs: type, length(if string), num_values
   if( (node = get_odb_node(path)) == NULL ){
      printf("get_odb_array:path[%s] not found\n", path); return(-1);
   }
   if( (attr = get_attr(node->attr, "type")) == NULL ){
      printf("get_odb_array:type attr not found\n", path); return(-1);
   }
   if( (tmp = attr_value(attr)) == NULL ){
      printf("get_odb_array:non-pointer-type-attr[%d]\n", attr->type);
      return(-1);
   }
   if( (*type = arraystr_idx(tmp, (char *)odb_typestr, DTYPE_STRLEN, ODB_DATATYPES)) == -1 ){ printf("get_odb_array:unknown datatype[%s]\n", tmp);
      return(-1);
   }
   if( (attr = get_attr(node->attr, "num_values")) == NULL ){
      printf("get_odb_array:num_values attr not found\n", path); return(-1);
   }
   if( (tmp = attr_value(attr)) == NULL ){
      printf("get_odb_array:non-pointer-num_values-attr[%d]\n", attr->type);
      return(-1);
   }
   if( sscanf(tmp, "%d", nvalues) < 1 ){
      printf("get_odb_array:cant read num_values[%s]\n", tmp);
      return(-1);
   }
   if( *type == 11 ){ // string
      if( (attr = get_attr(node->attr, "size")) == NULL ){
         printf("get_odb_array:string length attr not found\n", path);
         return(-1);
      }
      if( (tmp = attr_value(attr)) == NULL ){
         printf("get_odb_array:non-pointer-str-length-attr[%d]\n", attr->type);
         return(-1);
      }
      if( sscanf(tmp, "%d", size) < 1 ){
         printf("get_odb_array:cant read string-length[%s]\n", tmp);
         return(-1);
      }
      if( *size < 1 || *size >= MAX_ARRAYSTR_SIZE ){
         printf("get_odb_array:bad string-array length[%d]\n", *size);
         return(-1);
      }
   }
   if( (result = calloc(*nvalues, *size) ) == NULL ){
     printf("get_odb_array: failed malloc for %s\n", path); return(-1);
   }
   *value = i_array = f_array = result;
   if( (node = node->child_head) == NULL ){
      printf("get_odb_array:empty\n", path); return(-1);
   }
   while( node != NULL ){ tmp=node->name;
      if( strlen(tmp) != 5 || strncmp(tmp, "value", 5) != 0 ){
         printf("get_odb_array:non-value-key[%s]\n", tmp);
         node = node->next; continue;
      }
      if( (attr = get_attr(node->attr, "index")) == NULL ){
         printf("get_odb_array:index attr not found\n", path);
         node = node->next; continue;
      }
      if( (tmp = attr_value(attr)) == NULL ){
         printf("get_odb_array:non-ptr-index-attr-type[%d]\n", attr->type);
         node = node->next; continue;
      }
      if( sscanf(tmp, "%d", &index) < 1 ){
         printf("get_odb_array:cant read index[%s]\n", tmp);
         node = node->next; continue;
      }
      if( index < 0 || index >= *nvalues ){
         printf("get_odb_array:index[%d] out of range[0-%d]\n",index,*nvalues);
         node = node->next; continue;
      }
      if( (tmp = node_value(node)) == NULL ){
         printf("get_odb_array:non-pointer-value-type[%d]\n", attr->type);
         node = node->next; continue;
      }
      switch( *type ){
      case 0: case 1: case 2: case 3 :case 4: case 5: case 6: case 7: // ints
         if( sscanf(tmp, "%d", &i_val) < 1 ){
            printf("get_odb_array:cant read integer value[%s]\n", tmp);
            break;
         }
         ((int *)result)[index] = i_val; break;
      case 8: case 9: // floats
         if( sscanf(tmp, "%f", &f_val) < 1 ){
            printf("get_odb_array:cant read float value[%s]\n", tmp);
            break;
         }
         // doesnt work ?? ((float *)result)[index] = f_val; break;
         f_array[index] = f_val; break;
      case 11: // string
         if( strlen(tmp) >= *size ){
            printf("get_odb_array:string[%s] above specified length\n", tmp);
            break;
         }
         memcpy(result+index* *size, tmp, strlen(tmp) ); break;
      default:
         printf("get_odb_array:unhandled data type[%s]\n",odb_typestr[*type]);
         break;
      }
      node = node->next; continue;
   }
   return(0);
}

static Xml_val *get_odb_value(char *path, int *type)
{  // value is stored in value element, not subnode
   Xml_node *node = get_odb_node(path);
   if( node == NULL ){
      printf("get_odb_value:path[%s] not found\n", path); return(NULL);
   }
   *type = node->type;
   return(&node->value);
}

int read_odb_tree(int bank_len, int *bank_data)
{
   if( odb_tree != NULL ){ free_odb_tree(odb_tree); odb_tree=NULL; }
   read_xml_data((char *)bank_data, NULL, "root", 4, NULL, NULL);
   return(0);
}

// this is called recursively, when new start-tag is found
//    (recursion is to avoid having to maintain a list of parent nodes     )
//    (back to root or, even worse, having to search entire tree for parent)
char *read_xml_data(char *xml_file, Xml_node *current, char *name_ptr, int namelen, char *attr, char *attr_end)
{
   Xml_node *new_node, *parent = current;
   char *ptr = xml_file, *str, *ptr2;
   int len;
   
   if( (new_node = alloc_odb_element(name_ptr, namelen)) == NULL ){
      return(NULL);
   }
   if( current == NULL ){ // just starting - tree is currently empty
      current = odb_tree = new_node;
   } else { // all new elements are added to child-list
      if( current->child_head == NULL ){
         current = current->child_head = new_node;
      } else {
         current = current->child_head;
         while( current->next != NULL ){ current = current->next; }
         current = current->next = new_node;
      }
   }
   while( 1 ){ // add attrs [str points to end of tag]
      if( attr >= attr_end ){ break; }
      if( (str = strchr(attr,'=')) == NULL ){  // attr missing value - ignore this
         printf("read_xml_data: incomplete attr - ignoring attrs at %.16s\n", attr);
         break;
      }
      ptr2 = NULL;
      if( *(str+1) == '"'  ){ ptr2 = strchr(str+2,'"');  }
      if( *(str+1) == '\'' ){ ptr2 = strchr(str+2,'\''); }
      if( ptr2 == NULL ){  // attr missing value - ignore this
         printf("read_xml_data: incomplete attr - ignoring attrs at %.16s\n", attr);
         break;
      }
      insert_odb_attr(current, attr, str-attr, str+2, ptr2-str-2);
      attr = ptr2+2; // skip end-quote and space before next attr
   }
   while(1){ // have just started new key - read any data
      if( (str = strchr(ptr,'<')) == NULL ){ break; }
      if( str > ptr ){
         if( (current->value.p_val = malloc(str-ptr+1)) == NULL ){
            printf("read_xml_data: value malloc fail\n");
         }
         memcpy(current->value.p_val, ptr, str-ptr);
         current->value.p_val[str-ptr] = 0;
         current->type = XML_PTR;
         ptr = str;
      }
      // read next markup element
      if( !strncmp(ptr,"<?xml", 5 ) ){                // prolog - skip
         if( (str = strchr(ptr,'>')) == NULL ){ break; }
         while( isspace(*(++str)) ){ ; }
         ptr = str; continue;
      }
      if( !strncmp(ptr,"<!--", 4 ) ){                 // comment - skip
         if( (str = strchr(ptr,'>')) == NULL ){ break; }
         while( isspace(*(++str)) ){ ; }
         ptr = str; continue;
      }
      // start of new tag ...
      if( (str = strchr(ptr,'>')) == NULL ){
         printf("read_xml_data: incomplete tag at %.16s\n", ptr); return(NULL);
      }
      if( *(ptr+1) == '/' ){ // closing tag - check name
         len = str - ptr - 2;
         if( len != strlen(current->name) ||
             strncmp(ptr+2, current->name, len) != 0 ){
            printf("read_xml_data: mismatched start/end tags at %.16s\n", ptr); return(NULL);
         }
         return(str+1);
      }
      if( *(str-1)=='/' ){ ptr = str+1; continue; } // self-closing-tag - ignore
      if( (ptr2 = strchr(ptr,' ')) == NULL || ptr2 > str ){ // no attrs
         ptr = read_xml_data(str+1, current, ptr+1, str-ptr-1, NULL, NULL);
      } else {
         ptr = read_xml_data(str+1, current, ptr+1, ptr2-ptr-1, ptr2+1, str);
      }
      if( ptr == NULL ){ return(NULL); }
      continue;
   }
   return(ptr); // end of data (could also return null here?)
}

// ##########################################################################
// ##########################################################################

// find index of str in array of strings (-1 if not found)
int arraystr_idx(char *str, char *array, int size, int num)
{
   int i; char *ptr;
   for(i=0; i<num; i++){ ptr = array + size*i;
      if( strlen(str) == strlen(ptr) &&
          strncmp(str, ptr, strlen(str)) == 0 ){ break; }
   }
   return( (i<num) ? i : -1 );
}

int free_odb_attrlist(Xml_attr *current)
{
   Xml_attr *prev;
   while( current != NULL ){
      prev = current; current = current->next;
      if( prev->name != NULL ){ free(prev->name); } else {
         printf("free_odb_attrlist: NULL name pointer\n");
      }
      if( prev->type == XML_PTR ){
         if( prev->value.p_val != NULL ){ free(prev->value.p_val); } else {
            printf("free_odb_attrlist: NULL value pointer\n");
         }
      }
      free(prev);
   }
   return(0);
}
// strlen does not include trailing NULL
Xml_node *alloc_odb_element(char *name, int namelen)
{
   Xml_node *ptr;
   if( (ptr = malloc(sizeof(Xml_node))) == NULL){
      printf("alloc_odb_element: malloc fail\n"); return(NULL);
   }
   memset(ptr, 0, sizeof(Xml_node));
   if( (ptr->name = malloc( namelen+1 )) == NULL ){
      printf("alloc_odb_element: name malloc fail\n");
      free(ptr); ptr = NULL; return(NULL);
   }
   memcpy(ptr->name, name, namelen ); ptr->name[namelen] = 0;
   return(ptr);
}
int free_odb_element(Xml_node *current)
{
   if( current == NULL ){
      printf("free_odb_element: NULL pointer\n"); return(-1); }
   if( current->name != NULL ){ free(current->name ); } else {
      printf("free_odb_element: NULL name pointer\n");
   }
   if( current->type == XML_PTR ){
      if( current->value.p_val != NULL ){free(current->value.p_val); } else {
         printf("free_odb_element NULL value pointer\n");
      }  
   }
   if( current->attr != NULL ){ free_odb_attrlist(current->attr); }
   free(current);
   return(0);
}
// any nodes that have child_ptr are trees
int free_odb_tree(Xml_node *root)
{
   Xml_node *parent, *prev, *current;
   if( (current = root) == NULL ){ return(0); }
   while( 1 ){
      prev=current; current = current->next;
      if( prev->child_head != NULL ){
         free_odb_tree(prev->child_head);
      } else {
         free_odb_element(prev);
      }
      if( current == NULL ){ break; }
   }
   return(0);
}
int insert_odb_attr(Xml_node *current, char *name_ptr, int namelen, char *value, int vlen)
{
   Xml_attr *next, *attr = current->attr;
   if( (next = malloc(sizeof(Xml_attr))) == NULL ){
      printf("insert_odb_attr: attr malloc fail\n"); return(-1);
   }
   memset(next, 0, sizeof(Xml_attr));
   if( attr == NULL ){ attr = current->attr = next;
   } else {
      while( attr->next != NULL ){ attr = attr->next; }
      attr = attr->next = next;
   }
   if( (attr->name = malloc(namelen+1)) == NULL ){
      printf("insert_odb_attr: name malloc fail\n"); return(-1);
   }
   memcpy(attr->name, name_ptr, namelen);
   attr->name[namelen] = 0;
   if( (attr->value.p_val = malloc(vlen+1)) == NULL ){
      printf("insert_odb_attr: value malloc fail\n"); return(-1);
   }
   memcpy(attr->value.p_val, value, vlen);
   attr->value.p_val[vlen] = 0;
   attr->type = XML_PTR;
   return(0);
}
Xml_attr *get_attr(Xml_attr *attr, char *id)
{
   int len = strlen(id);
   while( attr != NULL ){
      if( strlen(attr->name) == len &&
          strncmp(attr->name, id, len) == 0 ){ break; }
      attr = attr->next;   
   }
   return(attr);
}
char *attr_value(Xml_attr *attr)
{
   char *tmp = attr->value.p_val;
   if( attr->type == XML_PTR && tmp != NULL ){ return(tmp); }
   return(NULL);
}
char *node_value(Xml_node *node)
{
   char *tmp = node->value.p_val;
   if( node->type == XML_PTR && tmp != NULL ){ return(tmp); }
   return(NULL);
}
// final path element can be dir, key or keyarray (others should all be dir)
Xml_node *get_odb_node(char *path)
{
   char *sep, *ptr = path, *tmp;
   Xml_node *node = odb_tree;
   Xml_attr *attr;
   int len;

   if( (node = odb_tree->child_head) == NULL){ // root
      printf("get_odb_node: empty odb tree\n"); return(NULL);
   }
   if( (node = node->child_head) == NULL){ // odb
      printf("get_odb_node: empty odb tree\n"); return(NULL);
   }
   while(1){
      if( (sep = strchr(ptr,'/')) == NULL ){ len = strlen(ptr);
      } else { len = sep-ptr; }
      while(1){ tmp = node->name;
         if( (strlen(tmp) == 3 && (strncmp(tmp, "dir", 3) == 0 ||
                                   strncmp(tmp, "key", 3) == 0) )   ||
             (strlen(tmp) == 8 && strncmp(tmp, "keyarray", 8) == 0) ){
            if( tmp[0] == 'k' && sep != NULL ){
               printf("get_odb_node: key mid path\n");
            }
            if( (attr = get_attr(node->attr, "name")) != NULL ){
               if( attr->type = XML_PTR && attr->value.p_val != NULL ){
                  if( strlen(attr->value.p_val) == len &&
                     strncmp(attr->value.p_val, ptr, len) == 0 ){
                     break;
                  }
               }
            }
         }
         if( (node = node->next) == NULL ){
            printf("get_odb_node: %s not found at %s\n", path, ptr); return(NULL);
         }
      }
      if( sep == NULL || strlen(sep+1) == 0 ){ break; }
      if( (node = node->child_head) == NULL){
         printf("get_odb_node: [%s] not found at %s\n", path, ptr); return(NULL);
      }
      ptr = sep+1;
   }
   return(node);
}
