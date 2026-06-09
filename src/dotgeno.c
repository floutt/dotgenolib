#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <stdbool.h>
#include <math.h>
#include "khash.h"
#include "dotgeno.h"

// macros
#define MAX(x, y) (((x) > (y)) ? (x) : (y))

#define SNP_NUM_COLS 6
#define SNP_VAR_COL 0
#define SNP_CHR_COL 1
#define SNP_CM_COL 2
#define SNP_POS_COL 3
#define SNP_REF_COL 4
#define SNP_ALT_COL 5

#define IND_NUM_COLS 3
#define IND_ID_COL 0
#define IND_SEX_COL 1
#define IND_POP_COL 2

#define BITS_IN_BYTE 8
#define RECORD_ELEM_SIZE_BITS 2
#define RECORD_ELEMS_PER_BYTE 4
#define RECORD_ELEMS_MASK_BASE 3

/* BASIC FUNCTIONS */
FILE* safe_read(char* filename, char* mode) {
	FILE* fp = fopen(filename, mode);
	if (fp == NULL) {
		fprintf(stderr, "ERROR: cannot open file %s\n", filename);
		exit(EXIT_FAILURE);
	}
	return fp;
}

size_t get_filesize(char* filename) {
	FILE* fp = safe_read(filename, "rb");
	fseek(fp, 0, SEEK_END);
	size_t file_size = ftell(fp);
	fclose(fp);
	return file_size;
}

uint64_t get_number_of_lines(char* filename) {
	uint64_t num_lines = 0;
	FILE* fp = safe_read(filename, "r");
	while(!feof(fp)) {
		num_lines += (fgetc(fp) == '\n');
	}
	fclose(fp);
	return num_lines;
}

char** get_column_elems(char* col_str, size_t ncol) {
	char* cur_col = strtok(col_str, " \n\t\r");
	size_t col_num = 0;
	char** col_elems = (char**) malloc(ncol * sizeof(char*));
	while(cur_col != NULL) {
		if(col_num == ncol) {
			if(strcmp(cur_col, "") == 0) {
				return col_elems;
			} else {
				fprintf(stderr, "ERROR: Too many columns in '.snp' or '.ind' file.\n");
				exit(EXIT_FAILURE);
			}
		} else {
			col_elems[col_num] = cur_col;
			col_num += 1;
			cur_col = strtok(NULL, " \n\t\r");
		}
	}
	if (col_num == ncol) {
		return col_elems;
	} else {
		fprintf(stderr, "ERROR: Too few columns in '.snp' or '.ind' file.\n");
		exit(EXIT_FAILURE);
	}
}

snp_data read_snp_file(char* filename) {
	FILE* fp = safe_read(filename, "r");
	uint64_t num_snps = get_number_of_lines(filename);
	char* line = NULL;
	size_t size = 0;
	ssize_t nread;
	
	snp_data snp_info = {
		.length = num_snps,
		.var_id = (char**) malloc(num_snps * sizeof(char*)),
		.chr = (char**) malloc(num_snps * sizeof(char*)),
		.pos_morgans = (double*) malloc(num_snps * sizeof(double)),
		.pos = (uint64_t*) malloc(num_snps * sizeof(uint64_t)),
		.ref = (char**) malloc(num_snps * sizeof(char*)),
		.alt = (char**) malloc(num_snps * sizeof(char*)),
		.rev_idx = kh_init(ID_MAP_STR),
		.hash = 0
	};

	size_t idx = 0;
	while((nread = getline(&line, &size, fp)) != -1) {
		char** elems = get_column_elems(line, SNP_NUM_COLS);
		snp_info.var_id[idx] = strdup(elems[SNP_VAR_COL]);
		snp_info.chr[idx] = strdup(elems[SNP_CHR_COL]);
		sscanf(elems[SNP_CM_COL], "%lf", &snp_info.pos_morgans[idx]);
		snp_info.pos[idx] = atoi(elems[SNP_POS_COL]);
		snp_info.ref[idx] = strdup(elems[SNP_REF_COL]);
		snp_info.alt[idx] = strdup(elems[SNP_ALT_COL]);
		// calculate hash
		snp_info.hash *= 17;
		snp_info.hash = snp_info.hash ^ hash_str(snp_info.var_id[idx]);
		free(elems);
		int ret;
		khiter_t k = kh_put(ID_MAP_STR, snp_info.rev_idx, snp_info.var_id[idx], &ret);
		if(ret == -1) {
			fprintf(stderr, "Error: Failed to insert variant name into hash table!\n");
			exit(EXIT_FAILURE);
		} else if(ret == 0) {
			fprintf(stderr, "Error: Multiple instances of variant %s found in the .snp file.\n", snp_info.var_id[idx]);
			exit(EXIT_FAILURE);
		}
		kh_value(snp_info.rev_idx, k) = idx;
		idx++;
	}
	free(line);
	fclose(fp);
	return snp_info;
}

ind_data read_ind_file(char* filename) {
	FILE* fp = safe_read(filename, "r");
	uint64_t num_inds = get_number_of_lines(filename);
	char* line = NULL;
	size_t size = 0;
	ssize_t nread;
	
	ind_data ind_info = {
		.length = num_inds,
		.ind_id = (char**) malloc(num_inds * sizeof(char*)),
		.sex = (char**) malloc(num_inds * sizeof(char*)),
		.population = (char**) malloc(num_inds * sizeof(char*)),
		.rev_idx = kh_init(ID_MAP_IND),
		.hash = 0
	};

	size_t idx = 0;
	while((nread = getline(&line, &size, fp)) != -1) {
		char** elems = get_column_elems(line, IND_NUM_COLS);
		ind_info.ind_id[idx] = strdup(elems[IND_ID_COL]);
		ind_info.sex[idx] = strdup(elems[IND_SEX_COL]);
		ind_info.population[idx] = strdup(elems[IND_POP_COL]);
		free(elems);
		// calculate hash
		ind_info.hash *= 17;
		ind_info.hash = ind_info.hash ^ hash_str(ind_info.ind_id[idx]);
		int ret;
		ind_idx ind_idx_struct = (ind_idx){ind_info.ind_id[idx], ind_info.population[idx]};
		khiter_t k = kh_put(ID_MAP_IND, ind_info.rev_idx, ind_idx_struct, &ret);
		if(ret == -1) {
			fprintf(stderr, "Error: Failed to insert individual into hash table!\n");
			exit(EXIT_FAILURE);
		} else if(ret == 0) {
			fprintf(stderr, "Error: Multiple instances of individual (%s, %s) found in the .snp file.\n", ind_info.ind_id[idx], ind_info.population[idx]);
			exit(EXIT_FAILURE);
		}
		kh_value(ind_info.rev_idx, k) = idx;
		idx++;
	}
	free(line);
	fclose(fp);
	return ind_info;
}

short write_snp_data(snp_data* snp_info, char* filename) {
	FILE* fp = safe_read(filename, "w+");
	for(size_t i = 0; i < snp_info->length; i++) {
		fprintf(fp, "%s\t%s\t%lf\t%zu\t%s\t%s\n",
				snp_info->var_id[i], snp_info->chr[i], snp_info->pos_morgans[i], snp_info->pos[i], snp_info->ref[i], snp_info->alt[i]);
	}
	return fclose(fp);
}

short write_ind_data(ind_data* ind_info, char* filename) {
	FILE* fp = safe_read(filename, "w+");
	for(size_t i = 0; i < ind_info->length; i++) {
		fprintf(fp, "%s\t%s\t%s\n",
				ind_info->ind_id[i], ind_info->sex[i], ind_info->population[i]);
	}
	return fclose(fp);
}

short get_snp_idx(snp_data* snp_info, char* var_name, size_t* idx) {
	khint_t k = kh_get(ID_MAP_STR, snp_info->rev_idx, var_name);
	if(k == kh_end(snp_info->rev_idx)) {
		return -1;
	} else {
		*idx = kh_value(snp_info->rev_idx, k);
		return 0;
	}
}

short get_ind_idx(ind_data* ind_info, char* ind_id, char* ind_pop, size_t* idx) {
	ind_idx iidx;
	iidx.ind_id = ind_id;
	iidx.ind_pop = ind_pop; 
	khint_t k = kh_get(ID_MAP_IND, ind_info->rev_idx, iidx);
	if(k == kh_end(ind_info->rev_idx)) {
		return -1;
	} else {
		*idx = kh_value(ind_info->rev_idx, k);
		return 0;
	}
}

pam_file_reader open_pam(char* filename, snp_data* snp_info, ind_data* ind_info) {
	pam_file_reader pf;
	size_t file_size = get_filesize(filename);
	pf.is_hdr_read = false;
	pf.is_open = true;
	pf.idx = 0;
	pf.n_ind = ind_info->length;
	pf.n_snp = snp_info->length;
	pf.record_size = file_size / (snp_info->length + 1);
	if ((pf.record_size * (snp_info->length + 1)) != file_size) {
		fprintf(stderr, "Invalid PACKEDANCESTRYMAP file. File size must be a multiple of the number of SNPs.\n");
		exit(EXIT_FAILURE);
	}
	pf.fp = safe_read(filename, "rb");
	return pf;
}

egn_file_reader open_egn(char* filename, snp_data* snp_info, ind_data* ind_info) {
	egn_file_reader ef;
	ef.is_open = true;
	ef.n_ind = ind_info->length;
	ef.n_snp = snp_info->length;
	size_t file_size = get_filesize(filename);
	ef.idx = 0;
	if (((ind_info->length + 1) * snp_info->length) != file_size) {
		fprintf(stderr, "ERROR: Invalid EIGENSTRAT file. Expected file of size %u bytes.\n", (ind_info->length + 1) * snp_info->length);
		exit(EXIT_FAILURE);
	}
	ef.fp = safe_read(filename, "rb");
	return ef;
}

int close_pam_file_reader(pam_file_reader* pf) {
	pf->is_open = false;
	return fclose(pf->fp);
}

int close_egn_file_reader(egn_file_reader* ef) {
	ef->is_open = false;
	return fclose(ef->fp);
}

hdr_data read_pam_header(pam_file_reader* pf) {
	if(pf->is_hdr_read) {
		fprintf(stderr, "Header already read!");
		exit(EXIT_FAILURE);
	}
	if(!pf->is_open) {
		fprintf(stderr, "ERROR: pam_file_reader has been closed.\n");
		exit(EXIT_FAILURE);
	}
	hdr_data hdr_info;
	char* hdr = (char*) malloc(sizeof(char) * pf->record_size);
	fread(hdr, 1, pf->record_size, pf->fp);
	sscanf(hdr, "GENO   %u %u %x %x", &hdr_info.n_ind, &hdr_info.n_snp, &hdr_info.ind_hash, &hdr_info.snp_hash);
	free(hdr);
	pf->is_hdr_read = true;
	return hdr_info;
}

uint8_t* read_pam_record(pam_file_reader* pf) {
	if(!pf->is_hdr_read) {
		fprintf(stderr, "Header must be read before reading records!\n");
		exit(EXIT_FAILURE);
	}
	if(!pf->is_open) {
		fprintf(stderr, "ERROR: pam_file_reader has been closed.\n");
		exit(EXIT_FAILURE);
	}
	if(pf->idx == pf->n_snp) {
		return NULL;
	}
	uint8_t* record = (uint8_t*)malloc(pf->n_ind * sizeof(uint8_t));
	size_t num_leftover_bytes = pf->record_size - (int)ceil((float)(pf->n_ind*RECORD_ELEM_SIZE_BITS) / BITS_IN_BYTE);
	uint8_t record_byte;
	for(int i = 0; i < pf->n_ind; i++) {
		uint8_t elem_pos = i%RECORD_ELEMS_PER_BYTE;
		uint8_t shift_by = (BITS_IN_BYTE-RECORD_ELEM_SIZE_BITS) - (RECORD_ELEM_SIZE_BITS*elem_pos);
		if(elem_pos == 0) {
			record_byte = getc(pf->fp);
		}
		record[i] = (record_byte & (RECORD_ELEMS_MASK_BASE << shift_by)) >> shift_by;
	}
	pf->idx += 1;
	fseek(pf->fp, num_leftover_bytes, SEEK_CUR);  // skip useless bytes
	return record;
}

uint8_t* read_egn_record(egn_file_reader* ef) {
	if(!ef->is_open) {
		fprintf(stderr, "ERROR: egn_file_reader has been closed.\n");
		exit(EXIT_FAILURE);
	}
	if(ef->idx == ef->n_snp) {
		return NULL;
	}
	uint8_t* record = (uint8_t*)malloc(ef->n_ind * sizeof(uint8_t));
	char* line = NULL;
	size_t size = 0;
	ssize_t nread = getline(&line, &size, ef->fp);
	if (nread != (ef->n_ind+1)) {
		fprintf(stderr, "ERROR: improperly formatted EIGENSTRAT geno file. Expected %u entries in line, got % u.\n", ef->n_ind+1, nread);
		exit(EXIT_FAILURE);
	}
	for(int i = 0; i < ef->n_ind; i++) {
		// check if valid character
		if(!((line[i] == '0') || (line[i] == '1') || (line[i] == '2') || (line[i] == '9'))) {
			fprintf(stderr, "ERROR: characters must be '0', '1', '2', or '9'.\n");
			exit(EXIT_FAILURE);
		}
		record[i] = ((line[i] == '0') * 0) + ((line[i] == '1') * 1) + ((line[i] == '2') * 2) + ((line[i] == '9') * NAN_VAL);
	}
	free(line);
	ef->idx += 1;
	return record;
}

void goto_var_egn(egn_file_reader* ef, snp_data* snp_info, char* var_name) {
	if(!ef->is_open) {
		fprintf(stderr, "ERROR: egn_file_reader closed.\n");
		exit(EXIT_FAILURE);
	}
	khint_t k = kh_get(ID_MAP_STR, snp_info->rev_idx, var_name);
	if(k == kh_end(snp_info->rev_idx)) {
		fprintf(stderr, "ERROR: variant %s is not in the .snp file!\n", var_name);
		exit(EXIT_FAILURE);
	}
	size_t idx_go = kh_value(snp_info->rev_idx, k);
	fseek(ef->fp, idx_go * (ef->n_ind+1), SEEK_SET);
	ef->idx = idx_go;
}

void goto_var_pam(pam_file_reader* pf, snp_data* snp_info, char* var_name) {
	if(!pf->is_open) {
		fprintf(stderr, "ERROR: pam_file_reader closed.\n");
		exit(EXIT_FAILURE);
	}
	if(!pf->is_hdr_read) {
		fprintf(stderr, "ERROR: header must be read before going to variant.\n");
		exit(EXIT_FAILURE);	
	}
	khint_t k = kh_get(ID_MAP_STR, snp_info->rev_idx, var_name);
	if(k == kh_end(snp_info->rev_idx)) {
		fprintf(stderr, "ERROR: variant %s is not in the .snp file!\n", var_name);
		exit(EXIT_FAILURE);
	}
	size_t idx_go = kh_value(snp_info->rev_idx, k);
	fseek(pf->fp, (idx_go + 1) * (pf->record_size), SEEK_SET);
	pf->idx = idx_go;
}

pam_file_writer pam_file_writer_init(char* filename, snp_data* snp_info, ind_data* ind_info) {
	pam_file_writer pfw;
	pfw.hdr_read = false;
	pfw.is_open = true;
	pfw.n_snp = snp_info->length;
	pfw.n_ind = ind_info->length;
	pfw.n_written_snp = 0;
	pfw.fp = safe_read(filename, "wb+");
	return pfw;
}

egn_file_writer egn_file_writer_init(char* filename, snp_data* snp_info, ind_data* ind_info) {
	egn_file_writer efw;
	efw.is_open = true;
	efw.n_snp = snp_info->length;
	efw.n_ind = ind_info->length;
	efw.n_written_snp = 0;
	efw.fp = safe_read(filename, "w+");
	return efw;
}

void write_pam_header(pam_file_writer* pfw, snp_data* snp_info, ind_data* ind_info) {
	if(pfw->hdr_read) {
		fprintf(stderr, "ERROR: pam_file_writer object header has already been written!\n");
		exit(EXIT_FAILURE);
	}
	int num_chars = fprintf(pfw->fp, "GENO   %u %u %x %x", 	pfw->n_ind, pfw->n_snp, ind_info->hash, snp_info->hash);
	if(num_chars < 0) {
		fprintf(stderr, "ERROR: unsuccessful write to PACKEDANCESTRYMAP file.\n");
		exit(EXIT_FAILURE);
	}
	size_t record_size = MAX(num_chars, (int)ceil((float)(pfw->n_ind * RECORD_ELEM_SIZE_BITS) / BITS_IN_BYTE));
  	size_t n_trailing_bytes_hdr	= record_size - num_chars;
	for(int i = 0; i < n_trailing_bytes_hdr; i++) {
		int ret = fputc('\0', pfw->fp);
		if(ret == EOF) {
			fprintf(stderr, "ERROR: unsuccessful write to PACKEDANCESTRYMAP file.\n");
			exit(EXIT_FAILURE);
		}
	}
	pfw->hdr_read = true;
	pfw->record_size = record_size;
}

void write_pam_record(pam_file_writer* pfw, uint8_t* dosages) {
	if(!pfw->hdr_read) {
		fprintf(stderr, "ERROR: header has not been written for PACKEDANCESTRYMAP file!\n");
		exit(EXIT_FAILURE);
	}
	if(!pfw->is_open) {
		fprintf(stderr, "ERROR: PACKEDANCESTRYMAP file is closed!\n");
		exit(EXIT_FAILURE);
	}
	if(pfw->n_written_snp >= pfw->n_snp) {
		fprintf(stderr, "ERROR: Cannot write to PACKEDANCESTRYMAP file. All SNP records have been written.\n");
		exit(EXIT_FAILURE);
	}
	uint8_t record_byte = 0;
	for(int i = 0; i < pfw->n_ind; i++) {
		if((i != 0) && ((i%(BITS_IN_BYTE/RECORD_ELEM_SIZE_BITS)) == 0)) {
			int ret = fputc(record_byte, pfw->fp);
			if(ret == EOF) {
				fprintf(stderr, "ERROR: unsuccessful write to PACKEDANCESTRYMAP file.\n");
				exit(EXIT_FAILURE);
			}
			record_byte = 0;
		}
		uint8_t elem_pos = i%RECORD_ELEMS_PER_BYTE;
		uint8_t shift_by = (BITS_IN_BYTE-RECORD_ELEM_SIZE_BITS) - (RECORD_ELEM_SIZE_BITS*elem_pos);
		record_byte = record_byte | ((RECORD_ELEMS_MASK_BASE << shift_by) & (dosages[i] << shift_by));
	}

	// write last record if not written
	if((pfw->n_ind % (BITS_IN_BYTE/RECORD_ELEM_SIZE_BITS)) != 0) {
		int ret = fputc(record_byte, pfw->fp);
		if(ret == EOF) {
			fprintf(stderr, "ERROR: unsuccessful write to PACKEDANCESTRYMAP file.\n");
			exit(EXIT_FAILURE);
		}
	}

	size_t n_trailing_bytes = pfw->record_size - (int)ceil((float)(pfw->n_ind * RECORD_ELEM_SIZE_BITS) / BITS_IN_BYTE);
	// write trailing bytes
	for(int i = 0; i < n_trailing_bytes; i++) {
		int ret = fputc('\0', pfw->fp);
		if(ret == EOF) {
			fprintf(stderr, "ERROR: unsuccessful write to PACKEDANCESTRYMAP file.\n");
			exit(EXIT_FAILURE);
		}
	}
	pfw->n_written_snp += 1;
}

void write_egn_record(egn_file_writer* efw, uint8_t* dosages) {
	if(!efw->is_open) {
		fprintf(stderr, "ERROR: EIGENSTRAT file is closed!\n");
		exit(EXIT_FAILURE);
	}
	if(efw->n_written_snp >= efw->n_snp) {
		fprintf(stderr, "ERROR: Cannot write to EIGENSTRAT file. All SNP records have been written.\n");
		exit(EXIT_FAILURE);
	}
	for(int i = 0; i < efw->n_ind; i++) {
		char dsg = ('0'*(dosages[i] == 0)) + ('1'*(dosages[i] == 1)) + ('2'*(dosages[i] == 2)) + ('9'*(dosages[i] == NAN_VAL));
		int ret = fputc(dsg, efw->fp);
		if(ret == EOF) {
			fprintf(stderr, "ERROR: unsuccessful write to EIGENSTRAT file.\n");
			exit(EXIT_FAILURE);
		}
	}
	int ret = fputc('\n', efw->fp);
	if(ret == EOF) {
		fprintf(stderr, "ERROR: unsuccessful write to EIGENSTRAT file.\n");
		exit(EXIT_FAILURE);
	}
	efw->n_written_snp += 1;
}

void free_str_array(char** arr, size_t length) {
	for(int i = 0; i < length; i++) {
		free(arr[i]);
	}
	free(arr);
}

void free_snp_data(snp_data* snp_info) {
	free_str_array(snp_info->var_id, snp_info->length);
	free_str_array(snp_info->chr, snp_info->length);
	free(snp_info->pos_morgans);
	free(snp_info->pos);
	free_str_array(snp_info->ref, snp_info->length);
	free_str_array(snp_info->alt, snp_info->length);
	kh_destroy(ID_MAP_STR, snp_info->rev_idx);
}

void free_ind_data(ind_data* ind_info) {
	free_str_array(ind_info->ind_id, ind_info->length);
	free_str_array(ind_info->sex, ind_info->length);
	free_str_array(ind_info->population, ind_info->length);
	kh_destroy(ID_MAP_IND, ind_info->rev_idx);
}
