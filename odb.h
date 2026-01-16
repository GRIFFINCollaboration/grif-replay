#define XML_INT    1
#define XML_SHORT  2
#define XML_FLOAT  3
#define XML_DOUBLE 4
#define XML_PTR    5
typedef union xml_value_union {
   int i_val;  short s_val;  char *p_val;
   float f_val; double d_val;
} Xml_val;
typedef struct xml_attr_struct Xml_attr;
struct xml_attr_struct {
   char *name; int type; Xml_val value; Xml_attr *next;
};
typedef struct xml_node_struct Xml_node;
// XML tree - top node contains rest of stuff => top->next should be NULL
//    array of same-parent-nodes stored as list using next pointer
//    => need second pointer to store any child-nodes-list
struct xml_node_struct {
   char *name; Xml_attr *attr; // attrs contain nxt pointer
   int type; Xml_val value;
   Xml_node *next;  Xml_node *child_head;
};
extern Xml_node *odb_tree;

int read_odb_tree(int bank_len, int *bank_data);
// basic data access ...
void *get_odb_value(char *path, int *type);
int get_odb_array(char *path, void **value, int *type, int *nvalues, int *size);
// more convenient data access with type conversion ...
int odbval_int(char *key, int *dst);
int odbval_float(char *key, float *dst);
int odbarray_int(char *key, int *dst, int maxlen);
int odbarray_float(char *key, float *dst, int maxlen);
