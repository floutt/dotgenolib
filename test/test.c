#include <stdlib.h>
#include <time.h>
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

// not needed when using generate_test_runner.rb
int main(void) {
    UNITY_BEGIN();
    RUN_TEST(test_equality_egn_pam);
	RUN_TEST(test_goto_equality_pam_egn);
    return UNITY_END();
}
