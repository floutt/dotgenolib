#define  _POSIX_C_SOURCE 200809L
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <sys/queue.h>
#include "unity.h"
#include "dotgeno.h"


void setUp(void) {
}

void tearDown(void) {
}

size_t get_random_index(size_t length) {
	srand(time(NULL));
	return((size_t)rand() % length);
}

void test_equality_egn_pam(void) {
	char* file_base_egn = "test/data/PAM_EGN/egn_example_%i.%s";
	char* file_base_pam = "test/data/PAM_EGN/pam_example_%i.%s";
	char buf[38];
    int n_sets = 20;
	for(int i = 1; i <= n_sets; i++) {
		sprintf(buf, file_base_pam, i, "snp");
		snp_data snp_info_pam = read_snp_file(buf);
		sprintf(buf, file_base_pam, i, "ind");
		ind_data ind_info_pam = read_ind_file(buf);
		sprintf(buf, file_base_pam, i, "geno");
		pam_file_reader pfr = pam_file_reader_init(buf, &snp_info_pam, &ind_info_pam);
		read_pam_header(&pfr);

		sprintf(buf, file_base_egn, i, "snp");
		snp_data snp_info_egn = read_snp_file(buf);
		sprintf(buf, file_base_egn, i, "ind");
		ind_data ind_info_egn = read_ind_file(buf);
		sprintf(buf, file_base_egn, i, "geno");
		egn_file_reader efr = egn_file_reader_init(buf, &snp_info_egn, &ind_info_egn);

		TEST_ASSERT_EQUAL_UINT(snp_info_pam.length, snp_info_egn.length);
		TEST_ASSERT_EQUAL_UINT(ind_info_pam.length, ind_info_egn.length);
		
		uint8_t* pam_record;
		uint8_t* egn_record;
		while(1) {
			pam_record = read_pam_record(&pfr);
			egn_record = read_egn_record(&efr);
			if(pam_record == NULL) {
				if(egn_record) {
					free(egn_record);
					TEST_FAIL();
				}
				break;
			}
			TEST_ASSERT_EQUAL_UINT8_ARRAY(pam_record, egn_record, ind_info_pam.length);
			free(pam_record);
			free(egn_record);
		}
		free_snp_data(&snp_info_pam);
		free_snp_data(&snp_info_egn);
		free_ind_data(&ind_info_pam);
		free_ind_data(&ind_info_egn);
		close_pam_file_reader(&pfr);
		close_egn_file_reader(&efr);
	}
}

void test_goto_equality_pam_egn(void) {
	char* file_base_egn = "test/data/PAM_EGN/egn_example_%i.%s";
	char* file_base_pam = "test/data/PAM_EGN/pam_example_%i.%s";
	char buf[38];
    int n_sets = 20;
	for(int i = 1; i <= n_sets; i++) {
		sprintf(buf, file_base_pam, i, "snp");
		snp_data snp_info_pam = read_snp_file(buf);
		sprintf(buf, file_base_pam, i, "ind");
		ind_data ind_info_pam = read_ind_file(buf);
		sprintf(buf, file_base_pam, i, "geno");
		pam_file_reader pfr = pam_file_reader_init(buf, &snp_info_pam, &ind_info_pam);
		read_pam_header(&pfr);

		sprintf(buf, file_base_egn, i, "snp");
		snp_data snp_info_egn = read_snp_file(buf);
		sprintf(buf, file_base_egn, i, "ind");
		ind_data ind_info_egn = read_ind_file(buf);
		sprintf(buf, file_base_egn, i, "geno");
		egn_file_reader efr = egn_file_reader_init(buf, &snp_info_egn, &ind_info_egn);

		TEST_ASSERT_EQUAL_UINT(snp_info_pam.length, snp_info_egn.length);
		TEST_ASSERT_EQUAL_UINT(ind_info_pam.length, ind_info_egn.length);
		
		int n_iters = 7;
		uint8_t* pam_record;
		uint8_t* egn_record;
		for(int j = 0; j < n_iters; j++) {
			short ret;
			ret = goto_var_pam(&pfr, &snp_info_pam, snp_info_pam.var_id[get_random_index(snp_info_pam.length)]);
			if(ret == -1) { TEST_FAIL(); };
			ret = goto_var_egn(&efr, &snp_info_egn, snp_info_egn.var_id[get_random_index(snp_info_egn.length)]);
			if(ret == -1) { TEST_FAIL(); };
			pam_record = read_pam_record(&pfr);
			egn_record = read_egn_record(&efr);
			TEST_ASSERT_EQUAL_UINT8_ARRAY(pam_record, egn_record, ind_info_pam.length);
			free(pam_record);
			free(egn_record);
		}
		free_snp_data(&snp_info_pam);
		free_snp_data(&snp_info_egn);
		free_ind_data(&ind_info_pam);
		free_ind_data(&ind_info_egn);
		close_pam_file_reader(&pfr);
		close_egn_file_reader(&efr);
	}
}

void test_ind_hash(void) {
	char* buf = NULL;
	size_t bufsize = 0;
	FILE* fp = fopen("test/data/hash_test/total/ind", "r");
	while(getline(&buf, &bufsize, fp) != -1) {
		char ind_file[37];
		uint32_t exp_hash;
		sscanf(buf, "%s %x\n", ind_file, &exp_hash);
		ind_data ind_info = read_ind_file(ind_file);
		TEST_ASSERT_EQUAL_UINT32(exp_hash, ind_info.hash);
		free_ind_data(&ind_info);
	}
	free(buf);
}

void test_snp_hash(void) {
	char* buf = NULL;
	size_t bufsize = 0;
	FILE* fp = fopen("test/data/hash_test/total/snp", "r");
	while(getline(&buf, &bufsize, fp) != -1) {
		char snp_file[37];
		uint32_t exp_hash;
		sscanf(buf, "%s %x\n", snp_file, &exp_hash);
		snp_data snp_info = read_snp_file(snp_file);
		TEST_ASSERT_EQUAL_UINT32(exp_hash, snp_info.hash);
		free_snp_data(&snp_info);
	}
	free(buf);
}

void test_ind_filtered_hash(void) {
	char* buf = NULL;
	size_t bufsize = 0;
	FILE* fp = fopen("test/data/hash_test/filtered/ind", "r");
	while(getline(&buf, &bufsize, fp) != -1) {
		char ind_file[37];
		uint32_t exp_hash;
		char* idx_str = malloc(sizeof(char) * ((bufsize - 46) + 1));
		sscanf(buf, "%s\t%s\t%x\n", ind_file, idx_str, &exp_hash);
		ind_data ind_info = read_ind_file(ind_file);
		
		// init idx list
		struct idx_head head;
		STAILQ_INIT(&head);
		char* elem = strtok(idx_str, ","); 
		// add elems
		while(elem) {
			struct idx_node* node = malloc(sizeof(struct idx_node));
			node->idx = atoi(elem);
			STAILQ_INSERT_TAIL(&head, node, nodes);
			elem = strtok(NULL, ",");
		}
		ind_data ind_filt;
		short ret = filter_ind_data(&ind_info, &ind_filt, &head);
		if(ret == -1) { TEST_FAIL(); }
		TEST_ASSERT_EQUAL_UINT32(exp_hash, ind_filt.hash);
		free_ind_data(&ind_info);
		free(idx_str);
		free_idx_list(&head);
	}
	free(buf);
}

void test_snp_filtered_hash(void) {
	char* buf = NULL;
	size_t bufsize = 0;
	FILE* fp = fopen("test/data/hash_test/filtered/snp", "r");
	while(getline(&buf, &bufsize, fp) != -1) {
		char snp_file[37];
		uint32_t exp_hash;
		char* idx_str = malloc(sizeof(char) * ((bufsize - 46) + 1));
		sscanf(buf, "%s\t%s\t%x\n", snp_file, idx_str, &exp_hash);
		snp_data snp_info = read_snp_file(snp_file);
		
		// init idx list
		struct idx_head head;
		STAILQ_INIT(&head);
		char* elem = strtok(idx_str, ","); 
		// add elems
		while(elem) {
			struct idx_node* node = malloc(sizeof(struct idx_node));
			node->idx = atoi(elem);
			STAILQ_INSERT_TAIL(&head, node, nodes);
			elem = strtok(NULL, ",");
		}
		snp_data snp_filt;
		short ret = filter_snp_data(&snp_info, &snp_filt, &head);
		if(ret == -1) { TEST_FAIL(); }
		TEST_ASSERT_EQUAL_UINT32(exp_hash, snp_filt.hash);
		free_snp_data(&snp_info);
		free(idx_str);
		free_idx_list(&head);
	}
	free(buf);
}

// not needed when using generate_test_runner.rb
int main(void) {
    UNITY_BEGIN();
    RUN_TEST(test_equality_egn_pam);
	RUN_TEST(test_goto_equality_pam_egn);
	RUN_TEST(test_ind_hash);
	RUN_TEST(test_snp_hash);
	RUN_TEST(test_ind_filtered_hash);
	RUN_TEST(test_snp_filtered_hash);
    return UNITY_END();
}
