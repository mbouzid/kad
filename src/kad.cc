#include <zlib.h>
#include <iostream>
#include <rocksdb/db.h>
#include <rocksdb/slice.h>
#include <rocksdb/write_batch.h>
#include <rocksdb/comparator.h>
#include <cassert>
#include <stdlib.h>
#include <math.h> // floor()
#include <cinttypes>
#include <sys/stat.h> // mkdir()
#include <sys/param.h> // MAXPATHLEN

#include "kseq.h"
#include "kstring.h"

KSEQ_INIT(gzFile, gzread)

#define KAD_VERSION "0.0.4"
#define KAD_DB_PREFIX ".kad"

#define KMER_LENGTH 32
#define NB_KMERS_PRINT 1000000
#define BUFFER_SIZE 10000

enum DNA_MAP {A, C, G, T};  // A=1, C=0, T=2, G=3
static const char NUCLEOTIDES[4] = { 'A', 'C', 'G', 'T' };

struct stat sb; // Use to check if files/directories exists

using namespace std;

uint64_t str_to_int(char* str)
{
  uint64_t strint = 0;
  for (size_t i = 0; i < KMER_LENGTH; i++) {
    uint8_t curr = 0;
    switch (str[i]) {
      case 'A': { curr = DNA_MAP::A; break; }
      case 'T': { curr = DNA_MAP::T; break; }
      case 'C': { curr = DNA_MAP::C; break; }
      case 'G': { curr = DNA_MAP::G; break; }
    }
    strint = strint << 2;
    strint = strint | curr;
  }
  return strint;
}

string int_to_str(uint64_t kmer)
{
  uint8_t base;
  string str;
  for (int i=KMER_LENGTH; i>0; i--) {
    base = (kmer >> (i*2-2)) & 3ULL;
    char chr;
		switch(base) {
			case DNA_MAP::A: { chr = 'A'; break; }
			case DNA_MAP::T: { chr = 'T'; break; }
			case DNA_MAP::C: { chr = 'C'; break; }
			case DNA_MAP::G: { chr = 'G'; break; }
		}
    str.push_back(chr);
  }
  return str;
}

typedef struct {
  uint16_t id;
  uint16_t n;
} count_t;

typedef struct {
  rocksdb::DB* samples_db;
  rocksdb::DB* counts_db;
} kad_db_t;

class KmerKeyComparator : public rocksdb::Comparator {
  public:
    int Compare(const rocksdb::Slice& a, const rocksdb::Slice& b) const {
      uint64_t kmer_a;
      uint64_t kmer_b;
      memcpy(&kmer_a, a.data(), sizeof(uint64_t));
      memcpy(&kmer_b, b.data(), sizeof(uint64_t));
      if(kmer_a < kmer_b) {
        return -1;
      } else if(kmer_b < kmer_a) {
        return +1;
      }
      return 0;
    }
    const char* Name() const { return "KmerKeyComparator"; }
    void FindShortestSeparator(std::string*, const rocksdb::Slice&) const { }
    void FindShortSuccessor(std::string*) const { }
};



kad_db_t* kad_open(const char* path) {
  kad_db_t* kad_db = (kad_db_t*)malloc(sizeof(kad_db_t));
  rocksdb::Options options_counts;
  rocksdb::Options options_samples;
  options_counts.create_if_missing = true;
  options_samples.create_if_missing = true;

  // Set options for counts
  KmerKeyComparator *cmp_kmers = new KmerKeyComparator(); // FIXME This should be deleted
  options_counts.comparator = cmp_kmers;
  options_counts.max_open_files = 1000;

  string db_path = path;
  db_path += "/";
  db_path += KAD_DB_PREFIX;

  if (stat(db_path.c_str(), &sb) != 0){
    if (mkdir(db_path.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH) != 0) {
      cerr << "Failed to create KAD directory: " << db_path << endl;
      exit(1);
    } else {
      cerr << "Successfully created KAD directory: " << db_path << endl;
    }
  }
  else if (!S_ISDIR(sb.st_mode)) {
    exit(1);
  }
  
 rocksdb::Status status;
 
 status = rocksdb::DB::Open(options_samples, string(db_path + "/samples").c_str(), &kad_db->samples_db);
 if(!status.ok()) {
   cerr << "Failed to open samples database" << endl;
   exit(2);
 }

 status = rocksdb::DB::Open(options_counts, string(db_path + "/counts").c_str(), &kad_db->counts_db);
 if(!status.ok()) {
   cerr << "Failed to open counts database" << endl;
   exit(2);
 }

  return kad_db;
}

void kad_destroy(kad_db_t *db) {
  delete db->samples_db;
  delete db->counts_db;
  free(db);
}

uint16_t add_sample(kad_db_t* db, const char* sample_name){
  uint16_t nb_keys;
  string value;
  // FIXME Test that "sample_name" is different from _nb_keys
  // FIXME This should be atomic
  rocksdb::Status s = db->samples_db->Get(rocksdb::ReadOptions(), "_nb_keys", &value);
  if(s.ok()) {
    nb_keys = atoi(value.c_str());
  } else {
    nb_keys = 0;
  }
  rocksdb::Slice key((char*)&nb_keys, sizeof(uint16_t));
  s = db->samples_db->Put(rocksdb::WriteOptions(), key, sample_name);
  if(!s.ok()) {
    cerr << "failed to add sample to the database" << endl;
    exit(3);
  } 
  // Update the number of keys
  std::string nb_keys_str = std::to_string(nb_keys+1);
  s = db->samples_db->Put(rocksdb::WriteOptions(), "_nb_keys", nb_keys_str);
  if(!s.ok()) {
    cerr << "failed to update the number of keys in the samples database" << endl;
    exit(3);
  }
  return nb_keys;
}

string get_sample(kad_db_t* db, uint16_t id) {
  string value;
  rocksdb::Slice key((char*)&id, sizeof(uint16_t));
  rocksdb::Status s = db->samples_db->Get(rocksdb::ReadOptions(), key, &value);
  return value;
}

void print_counts(kad_db_t* db, size_t nb_counts, count_t* counts) {
  // FIXME add local hash for sample names
  for(size_t i = 0; i < nb_counts; i++) {
    string sample_name = get_sample(db, counts[i].id);
    if(i > 0)
      cout << "\t";
    cout << sample_name << "|" << counts[i].n;
  }
}

int kad_test(kad_db_t* db, int argc, char **argv) {
  //char kmer[33] = "AGAGGAGGGACGGGCTGAAAAAGTACTCATTG";
  return 0;
}

int kad_dump(kad_db_t* db, int argc, char **argv)
{
  rocksdb::Iterator* it = db->counts_db->NewIterator(rocksdb::ReadOptions());
  for (it->SeekToFirst(); it->Valid(); it->Next()) {
    uint64_t kmer_int_found;
    memcpy(&kmer_int_found, it->key().data(), sizeof(uint64_t));

    count_t *counts = (count_t*)it->value().data();
    size_t nb_counts = it->value().size() / sizeof(count_t);

    std::cout << int_to_str(kmer_int_found) << "\t";
    print_counts(db, nb_counts, counts);

    cout << endl;
  }
  return 0;
}

uint64_t rand_uint64(void) {
  uint64_t r = 0;
  for (int i=0; i<64; i += 30) {
    r = r*(RAND_MAX + (uint64_t)1) + rand();
  }
  return r;
}

int kad_random_query(kad_db_t* db, int argc, char **argv) {

  if (argc < 2) {
		fprintf(stderr, "\n");
		fprintf(stderr, "Usage:   kad random_query [options] nb_queries\n\n");
		return 1;
  }

  size_t nb_queries = atoi(argv[1]);
  string value;

  for(size_t i = 0; i < nb_queries; i++) {
    uint64_t random_kmer = rand_uint64();
    rocksdb::Slice key((char*)&random_kmer, sizeof(uint64_t));
    rocksdb::Status s = db->counts_db->Get(rocksdb::ReadOptions(), key, &value);
  }

  return 0;
}

int kad_query(kad_db_t* db, int argc, char **argv) {

  if (argc < 2) {
		fprintf(stderr, "\n");
		fprintf(stderr, "Usage:   kad query [options] kmer\n\n");
		return 1;
  }

  std::string value;
  char *kmer = argv[1];

  uint64_t kmer_int = str_to_int(kmer);

  rocksdb::Slice key((char*)&kmer_int, sizeof(uint64_t));

  rocksdb::Status s = db->counts_db->Get(rocksdb::ReadOptions(), key, &value);
  if (s.ok()) {
    cout << kmer << "\t";

    count_t *counts = (count_t*)value.data();
    size_t nb_counts = value.size() / sizeof(count_t);

    print_counts(db, nb_counts, counts);

    cout << endl;
  }

  return 0;
}

int kad_index(kad_db_t* db, int argc, char **argv)
{

  if (argc < 2) {
		fprintf(stderr, "\n");
		fprintf(stderr, "Usage:   kad index [options] sample_name counts.tsv\n\n");
		return 1;
  }

  std::string value;
  char *sample_name = argv[1];
  char *file = argv[2];

  uint16_t sample_id = add_sample(db, sample_name);

  gzFile fp;
	kstream_t *ks;
	kstring_t *str,*kmer;
  int dret;
  uint16_t count;
  uint64_t kmer_int;
  size_t nb_kmers = 0;
  size_t batch_i = 0;
  rocksdb::WriteBatch batch;

  kmer  = (kstring_t*)calloc(1, sizeof(kstring_t));
  str   = (kstring_t*)calloc(1, sizeof(kstring_t));
  fp = gzopen(file, "r");
  if(!fp) { fprintf(stderr, "Failed to open %s\n", file); exit(EXIT_FAILURE); }

  ks = ks_init(fp);

  while (ks_getuntil(ks, 0, str, &dret) >= 0) {
    kputs(str->s,kmer);
    if(dret != '\n') {
      if(ks_getuntil(ks, 0, str, &dret) > 0 && isdigit(str->s[0])) {
        count = (uint16_t)atoi(str->s);
        
        kmer_int = str_to_int(kmer->s);

        rocksdb::Slice key((char*)&kmer_int, sizeof(uint64_t));

        count_t *counts;
        size_t nb_counts = 1;

        rocksdb::Status s = db->counts_db->Get(rocksdb::ReadOptions(), key, &value);

        // The k-mer is already in the database
        if(s.ok()) {
          nb_counts = (value.size() / sizeof(count_t)) + 1;
          counts = (count_t*)malloc(sizeof(count_t) * nb_counts);
          memcpy(counts,value.data(),value.size());
        } else {
          counts = (count_t*)malloc(sizeof(count_t));
        }

        counts[nb_counts-1] = { sample_id, count };
        
        rocksdb::Slice counts_value((char*)counts, nb_counts * sizeof(count_t));

        if(batch_i < BUFFER_SIZE) {
          batch.Put(key, counts_value);
          ++batch_i;
        }

        if(batch_i == BUFFER_SIZE ) {
          s = db->counts_db->Write(rocksdb::WriteOptions(), &batch);
          if(!s.ok()) {
            cerr << s.ToString() << endl;
            exit(4);
          }
          batch.Clear();
          batch_i = 0;
        }
        free(counts);
      }
    }
    kmer->l = 0;
    nb_kmers++;
    if(nb_kmers % NB_KMERS_PRINT == 0)
      cerr << nb_kmers  << " kmers loaded" << endl;
  }

  cerr << "Successfully loaded " << nb_kmers << " kmers" << endl;

  ks_destroy(ks);
  gzclose(fp);
  free(str->s); free(str);
  free(kmer->s); free(kmer);
  return 0;
}

/* main function */
static int usage()
{
	fprintf(stderr, "\n");
	fprintf(stderr, "Usage:   kad <command> <arguments>\n");
	fprintf(stderr, "Version: %s\n\n", KAD_VERSION);
	fprintf(stderr, "Command: index      Index k-mer counts from a samples\n");
	fprintf(stderr, "         query      Query the KAD database\n");
	fprintf(stderr, "         dump       Dump the KAD database\n");
	fprintf(stderr, "\n");
	return 1;
}

int main(int argc, char *argv[])
{
	if (argc == 1) return usage();

  char cwd[MAXPATHLEN];
  getcwd(cwd, MAXPATHLEN);
  kad_db_t* db = kad_open(cwd);

	if (strcmp(argv[1], "index") == 0) kad_index(db, argc-1, argv+1);
  else if (strcmp(argv[1], "dump") == 0) kad_dump(db, argc-1, argv+1);
  else if (strcmp(argv[1], "query") == 0) kad_query(db, argc-1, argv+1);
  else if (strcmp(argv[1], "random_query") == 0) kad_random_query(db, argc-1, argv+1);
  else if (strcmp(argv[1], "test") == 0) kad_test(db, argc-1, argv+1);
	else {
		fprintf(stderr, "[main] unrecognized command '%s'. Abort!\n", argv[1]);
		return 1;
	}

  kad_destroy(db);
	return 0;
}
