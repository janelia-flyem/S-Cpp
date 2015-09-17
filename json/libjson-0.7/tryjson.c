#include "json.h"

#include <stdio.h>
int my_callback(void *userdata, int type, const char *data, uint32_t length)
{
	printf("Callback called\n");
	FILE *output = (userdata) ? (FILE *)userdata : stdout;
	switch (type) {
	case JSON_OBJECT_BEGIN:
	case JSON_ARRAY_BEGIN:
		fprintf(output, "entering %s\n", (type == JSON_ARRAY_BEGIN) ? "array" : "object");
		break;
	case JSON_OBJECT_END:
	case JSON_ARRAY_END:
		fprintf(output, "leaving %s\n", (type == JSON_ARRAY_END) ? "array" : "object");
		break;
	case JSON_KEY:
	case JSON_STRING:
	case JSON_INT:
	case JSON_FLOAT:
		fprintf(output, "value %*s\n", length, data);
		break;
	case JSON_NULL:
		fprintf(output, "constant null\n"); break;
	case JSON_TRUE:
		fprintf(output, "constant true\n"); break;
	case JSON_FALSE:
		fprintf(output, "constant false\n"); break;
	}
return 0;
}

int main(int argc, char **argv)
{

json_parser parser;
void *my_callback_data = NULL;

if (json_parser_init(&parser, NULL, my_callback, my_callback_data)) {
    fprintf(stderr, "something wrong happened during init\n");
    }

int len, ret;
char block[1024];

FILE *fd = fopen(argv[1],"r");
if (fd == NULL) {
    printf("Could not open\n");
    return 42;
    }
while ((len = fread(block, 1, 1024, fd)) > 0) {
    printf("len = %d\n", len);
    ret = json_parser_string(&parser, block, len, NULL);
    if (ret) {
	    /* error happened : print a message or something */
            printf("Syntax error\n");
	    break;
	}
    }
fclose(fd);

json_parser_free(&parser);
}
