#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <stdbool.h>

// macros
#define SNP_VAR_COL 0
#define SNP_CHR_COL 1
#define SNP_CM_COL 2
#define SNP_POS_COL 3
#define SNP_REF_COL 4
#define SNP_ALT_COL 5

typedef struct {
	size_t length;
	char** var_id;
	char** chr;
	double* cm;
	uint64_t* pos;
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

typedef struct {
	size_t start;
	size_t length;
} col_pos;

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
	if(curr == EOF) {
		return NULL;
	}
	size_t end_pos = ftell(fp);
	fseek(fp, start_pos, SEEK_SET);
	char* out_str = (char*)malloc(len+1);
	fgets(out_str, len+1, fp);
	fseek(fp, end_pos, SEEK_SET);
	return out_str;
}

uint64_t get_number_of_lines(char* filename) {
	uint64_t num_lines = 0;
	FILE *fp = fopen(filename, "r");
	while(!feof(fp)) {
		num_lines += (fgetc(fp) == '\n');
	}
	fclose(fp);
	return num_lines;
}

col_pos* get_column_lengths(char* col_str, uint64_t ncol) {
	bool on_whitespace = 1;
	size_t len = strlen(col_str);
	size_t num_cols = 0;
	size_t col_len = 1;
	col_pos* col_info = (col_pos*) malloc(ncol * sizeof(col_pos));
	for(int i = 0; i < len; i++) {
		if(num_cols > ncol) { perror("Too many columns"); }
		if (isspace(col_str[i])) {
			on_whitespace = 1;
			continue;
		} else if(on_whitespace) {
			on_whitespace = 0;
			if(num_cols > 0) {
				col_info[num_cols - 1].start = i - (col_len+1);
				col_info[num_cols - 1].length = col_len;
			}
			col_len = 1;
			num_cols += 1;
		} else {
			col_len += 1;
		}
	}
	if(num_cols < ncol) { perror("Too few columns"); }
	col_info[num_cols - 1].start =  len - col_len;
	col_info[num_cols - 1].length = col_len;
	return col_info;
}

snp_data read_snp_file(char* filename) {
	FILE *fp = fopen(filename, "r");
	uint64_t num_snps = get_number_of_lines(filename);
	snp_data snp_info = {
		.length = num_snps,
		.var_id = (char**) malloc(num_snps * sizeof(char*)),
		.chr = (char**) malloc(num_snps * sizeof(char*)),
		.cm = (double*) malloc(num_snps * sizeof(double)),
		.pos = (uint64_t*) malloc(num_snps * sizeof(uint64_t)),
		.ref = (char**) malloc(num_snps * sizeof(char*)),
		.alt = (char**) malloc(num_snps * sizeof(char*)),
		.hash = 0
	};

	//while(1) {

	//}
	return snp_info;
}
