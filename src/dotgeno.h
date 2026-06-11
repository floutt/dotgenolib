#include <string.h>
#include <sys/queue.h>
#include "khash.h"

/*! @brief NAN value encoding */
#define NAN_VAL 3


struct idx_node {
	size_t idx;
	STAILQ_ENTRY(idx_node) nodes;
};

STAILQ_HEAD(idx_head, idx_node);

typedef struct {
	char** strs;
	size_t length;
} string_array;

/*! @typedef
 * @brief Struct for indexing ind_data objects
 *
 * This is used to access information from ".ind" files.
 */
typedef struct {
	char* ind_id; /**< Intrapopulation identifier for individual */
	char* ind_pop; /**< Unique population identifier*/
} ind_idx;

// String linker for .ind hashing
char IND_LINK[20] = "gzvrEy55bcEN0gqRqvL6";
#define IND_LINK_LEN 20

/**
 *  @brief Hash function for ind_idx for use with khash.h, based on the djb2 hash function
 *
 *  @param[in] ind_idx_struct ind_idx structure to be hashed
 *
 *  @return hash value of type khint_t
 */
khint_t hash_ind_idx(ind_idx ind_idx_struct) {
	khint_t hash = 5381;
	
	// first string
	int len = strlen(ind_idx_struct.ind_id);
	for(int i = 0; i < len; i++) {
		hash = ((hash << 5) + hash) + ind_idx_struct.ind_id[i];
	}
	
	// linker string
	for(int i = 0; i < IND_LINK_LEN; i++) {
		hash = ((hash << 5) + hash) + IND_LINK[i];
	}

	// second string
	len = strlen(ind_idx_struct.ind_pop);
	for(int i = 0; i < len; i++) {
		hash = ((hash << 5) + hash) + ind_idx_struct.ind_pop[i];
	}
	return hash;
}

/**
 *  @brief Hash function for single string for use in evaluating integrity of .ind and .snp files
 *
 *  This function is part of the hashing process used for these files. This one operates on
 *  individual IDs or variant IDs. The hash is built as the .snp and .ind files are being read.
 *  See read_snp_file and read_ind_file for more details on their implementation.
 *
 *  @param[in] str string to be hashed
 *
 *  @return hash value of type uint32_t
 */
uint32_t hash_str(char* str) {
	uint32_t hash_out = 0;
	size_t str_len = strlen(str);
	for(int i = 0; i < str_len; i++) {
		hash_out *= 23;
		hash_out += str[i];
	}
	return(hash_out);
}

/*! @brief comparison function for ind_idx hash. For use with khash.h */
#define ind_idx_equal(a, b) ((strcmp((a).ind_id, (b).ind_id) == 0) && strcmp((a).ind_pop, (b).ind_pop) == 0)

/*! @brief string to size_t hash table for use with snp_data structs */
KHASH_MAP_INIT_STR(ID_MAP_STR, size_t)
/*! @brief ind_idx to size_t hash table for use with ind_data structs */ 
KHASH_INIT(ID_MAP_IND, ind_idx, size_t, true, hash_ind_idx, ind_idx_equal)

/*! @typedef
 * @brief struct containing information encoded in ".snp" files
 */
typedef struct {
	size_t length; /**< Total number of elements for each array in struct */
	char** var_id; /**< Array of variant identifiers */
	char** chr; /**< Array of chromosomes of variants */
	double* pos_morgans; /**< Array of genetic positions in Morgans */
	uint64_t* pos; /**< Array of genetic base pair positions */
	char** ref; /**< Array of reference alleles */
	char** alt; /**< Array of alternative alleles */ 
	khash_t(ID_MAP_STR)* rev_idx; /**< Hash table which returns index of the given var_id */
	uint32_t hash; /**< Hash value for struct */
} snp_data;

/*! @typedef
 * @brief struct containing information encoded in ".ind" files
 */
typedef struct {
	size_t length; /**< Total number of elements for each array in struct */
	char** ind_id; /**< Array of intrapopulation identifiers for individuals */
	char** sex; /**< Array of sex for individuals */
	char** population; /**< Array of sex for individuals */
	khash_t(ID_MAP_IND)* rev_idx; /**< Hash table which returns index of the given ind_idx struct */
	uint32_t hash; /**< Hash value for struct */
} ind_data;

/*! @typedef
 * @brief struct containing information necessary to read PACKEDANCESTRYMAP files
 */
typedef struct {
	bool is_hdr_read; /**< Whether the header of the PACKEDANCESTRYMAP file has been read */
	bool is_open; /**< Whether the header of the PACKEDANCESTRYMAP file has been opened */
	size_t idx; /**< Index of the current record being iterated */
	size_t record_size; /**< Size of each PACKEDANCESTRYMAP record */
	size_t n_ind; /**< Number of individuals in the file */
	size_t n_snp; /**< Number of genetic variants in the file */
	FILE* fp;  /**< File pointer to PACKEDANCESTRYMAP file */ 
} pam_file_reader;

/*! @typedef
 * @brief struct containing information necessary to read EIGENSTRAT files
 */
typedef struct {
	bool is_open; /**< Whether the header of the EIGENSTRAT file has been opened */
	size_t idx; /**< Index of the current record being iterated */
	size_t n_ind; /**< Number of individuals in the file */
	size_t n_snp; /**< Number of genetic variants in the file */  
	FILE* fp; /**< File pointer to EIGENSTRAT file */
} egn_file_reader;

typedef struct {
	bool hdr_read;
	bool is_open;
	size_t n_snp;
	size_t n_ind;
	size_t n_written_snp;
	size_t record_size;
	FILE* fp;
} pam_file_writer;

typedef struct {
	bool is_open;
	size_t n_snp;
	size_t n_ind;
	size_t n_written_snp;
	FILE* fp;
} egn_file_writer;

pam_file_writer pam_file_writer_init(char* filename, snp_data* snp_info, ind_data* ind_info);

egn_file_writer egn_file_writer_init(char* filename, snp_data* snp_info, ind_data* ind_info);

void write_pam_header(pam_file_writer* pfw, snp_data* snp_info, ind_data* ind_info);

void write_pam_record(pam_file_writer* pfw, uint8_t* dosages);

void write_egn_record(egn_file_writer* efw, uint8_t* dosages);

/*! @typedef
 * @brief struct containing information encoded in PACKEDANCESTRYMAP file headers
 */
typedef struct {
	size_t n_ind; /**< Number of individuals in the file */
	size_t n_snp; /**< Number of genetic variants in the file */ 
	uint32_t ind_hash; /**< Expected hash value for corresponding ".ind" file */
	uint32_t snp_hash; /**< Expected hash value for corresponding ".snp" file */
} hdr_data;

/**
 *  @brief reads information from .snp file and stores it in snp_data struct
 *
 *  @param[in] filename path to .snp file
 *
 *  @return snp_data struct containing information from .snp file
 */ 
snp_data read_snp_file(char* filename);

/**
 *  @brief reads information from .ind file and stores it in ind_data struct
 *
 *  @param[in] filename path to .ind file
 *
 *  @return ind_data struct containing information from .ind file
 */ 
ind_data read_ind_file(char* filename);

/**
 * @brief writes contents of snp_data object to file
 *
 * @param[in] snp_info pointer to snp_data object
 * @param[in] filename path to output file
 *
 * @return status code indicating whether or not file was successfully read and closed
 * @retval 0 not successfully closed 
 * @retval EOF file not successfully closed
 */
short write_snp_data(snp_data* snp_info, char* filename);

/**
 * @brief writes contents of ind_data object to file
 *
 * @param[in] ind_info pointer to ind_data object
 * @param[in] filename path to output file
 *
 * @return status code indicating whether or not file was successfully read and closed
 * @retval 0 not successfully closed 
 * @retval EOF file not successfully closed
 */
short write_ind_data(ind_data* ind_info, char* filename);

/**
 * @brief gets the index for snp with the identifier var_name
 *
 * @param[in] snp_info pointer to snp_data struct storing variant information
 * @param[in] var_name ID of desired variant
 * @param[out] idx pointer to where index will be stored
 *
 * @return status code indicating whether or not variant was found in snp_info
 * @retval -1 variant not found
 * @retval 0 variant found
 */
short get_snp_idx(snp_data* snp_info, char* var_name, size_t* idx);

/**
 * @brief gets the index for an individual with with certain individual ID AND population
 *
 * @param[in] ind_info pointer to ind_data struct storing variant information
 * @param[in] ind_id intrapopulation identifier for individual
 * @param[in] ind_pop population of individual
 * @param[out] idx pointer to where index will be stored
 *
 * @return status code indicating whether or not individual was found in ind_info
 * @retval -1 individual not found
 * @retval 0 individual found
 */
short get_ind_idx(ind_data* ind_info, char* ind_id, char* ind_pop, size_t* idx);

void get_multiple_snp_idx(snp_data* snp_info, string_array* var_names, struct idx_head* head);

void get_multiple_ind_idx(ind_data* ind_info, string_array* ind_ids, string_array* ind_pops, struct idx_head* head);

short filter_snp_data(snp_data* snp_in, snp_data* snp_out, struct idx_head* head);

short filter_ind_data(ind_data* ind_in, ind_data* ind_out, struct idx_head* head);

/**
 * @brief opens a PACKEDANCESTRYMAP file
 *
 * @param[in] filename path to PACKEDANCESTRYMAP file
 * @param[in] snp_info pointer to corresponding snp_data object
 * @param[in] ind_info pointer to corresponding ind_data object
 *
 * @return object with file pointer and information necessary to read PACKEDANCESTRYMAP file 
 */
pam_file_reader open_pam(char* filename, snp_data* snp_info, ind_data* ind_info);

/**
 * @brief opens a EIGENSTRAT file
 *
 * @param[in] filename path to EIGENSTRAT file
 * @param[in] snp_info pointer to corresponding snp_data object
 * @param[in] ind_info pointer to corresponding ind_data object
 *
 * @return object with file pointer and information necessary to read EIGENSTRAT file 
 */
egn_file_reader open_egn(char* filename, snp_data* snp_info, ind_data* ind_info);

/**
 * @brief closes pam_file_reader object
 *
 * @param[in,out] pf pointer to pam_file_reader object
 *
 * @return status code indicating whether or not the file was successfully closed
 * @retval 0 file successfully closed
 * @retval EOF file not successfully closed
 */
int close_pam_file_reader(pam_file_reader* pf);

/**
 * @brief closes egn_file_reader object
 *
 * @param[in,out] ef pointer to egn_file_reader object
 *
 * @return status code indicating whether or not the file was successfully closed
 * @retval 0 file successfully closed
 * @retval EOF file not successfully closed
 */
int close_egn_file_reader(egn_file_reader* ef);

/**
 * @brief reads header record of PACKEDANCESTRY map file.
 * This must be read before any other record is read
 *
 * @param[in,out] pf pointer to pam_file_reader object 
 *
 * @return hdr_data object containing important header information
 */
hdr_data read_pam_header(pam_file_reader* pf);

/**
 * @brief reads record of dosages of the current genetic variant from pam_file_reader object and iterates to the next one
 *
 * @param[in,out] pf pointer to pam_file_reader object
 *
 * @return array of allelic dosages of length pf->length. Dosages must be either 0,1,2, or NAN_VAL.
 */
uint8_t* read_pam_record(pam_file_reader* pf);

/**
 * @brief reads record of dosages of the current genetic variant from egn_file_reader object and iterates to the next one
 *
 * @param[in,out] ef pointer to egn_file_reader object
 *
 * @return array of allelic dosages of length ef->length. Dosages must be either 0,1,2, or NAN_VAL.
 */
uint8_t* read_egn_record(egn_file_reader* ef);

/**
 * @brief sends egn_file_reader object to the location of a particular variant
 *
 * @param[in,out] ef pointer to egn_file_reader object
 */
void goto_var_egn(egn_file_reader* ef, snp_data* snp_info, char* var_name);

/**
 * @brief sends pam_file_reader object to the location of a particular variant
 *
 * @param[in,out] pf pointer to pam_file object
 */
void goto_var_pam(pam_file_reader* pf, snp_data* snp_info, char* var_name);

/**
 * @brief frees snp_data object
 */
void free_snp_data(snp_data* snp_info);

/**
 * @brief frees ind_data object
 */
void free_ind_data(ind_data* ind_info);
