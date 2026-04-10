#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

typedef struct {
	size_t length;
	char** var_id;
	char** chr;
	double cm;
	uint64_t pos;
	char** ref;
	char** alt;
	uint32_t hash;
} snp_data;

typedef struct {
	size_t length;
	char** ind_id;
	char** sex;
	char** population;
	uint32_t hash;
} ind_data;

size_t get_filesize(char* filename) {
	FILE *fp = fopen(filename, "rb");
	fseek(fp, 0, SEEK_END);
	size_t file_size = ftell(fp);
	fclose(fp);
	return file_size;
}

char* read_line(FILE* fp) {
	/* get length of line, newline and EOF excluded */
	size_t start_pos = ftell(fp);
	char curr = fgetc(fp);
	uint64_t len = 0;
	while((curr != '\n') && (curr != EOF)) {
		len++;
		curr = fgetc(fp);
	}
	size_t end_pos = ftell(fp);
	fseek(fp, start_pos, SEEK_SET);
	char* out_str = (char*)malloc(len+1);
	fgets(out_str, len+1, fp);
	fseek(fp, end_pos, SEEK_SET);
	return out_str;
}

snp_data read_snp_file(char* filename) {
	FILE *fp = fopen(filename, "r");
}
