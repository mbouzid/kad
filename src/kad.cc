#include <zlib.h>
#include <iostream>
#include <rocksdb/db.h>
#include <rocksdb/slice.h>
#include <cassert>
#include <stdlib.h>
#include <jansson.h>
#include <math.h> // floor()
#include <cinttypes>

#include "kseq.h"
#include "kstring.h"

KSEQ_INIT(gzFile, gzread)

#define KAD_VERSION "0.0.3"
#define LEVELDB_PATH "/tmp/testdb"

#define KMER_LENGTH 32

enum DNA_MAP {C, A, T, G};  // A=1, C=0, T=2, G=3
static const char NUCLEOTIDES[4] = { 'C', 'A', 'T', 'G' };

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


int kad_test(int argc, char **argv) {
  char kmer[33] = "AGAGGAGGGACGGGCTGAAAAAGTACTCATTG";

  uint64_t kmer_int = str_to_int(kmer);

  string kmer_str = int_to_str(kmer_int);
  cout << kmer << "\t" << int_to_str(kmer_int) << endl;
  

  char *buf = (char*)malloc(sizeof(uint64_t));
  memcpy(buf, (char*)&kmer_int,sizeof(uint64_t));
  
  uint64_t kmer_int_found;
  memcpy(&kmer_int_found, buf, sizeof(uint64_t));

  cout << int_to_str(kmer_int) << endl;

  //char *buf = malloc(sizeof(kmer_int));
  //memcpy(buf, &kmer_int, sizeof(kmer_int));
  //Slice(buf, sizeof(kmer_int));

  //string key;
  //PutFixed64(&key,kmer_int);

  //char packed_kmer[KMER_PACKED_LENGTH];
  //char unpacked_kmer[KMER_LENGTH + 1];
  //unpacked_kmer[KMER_LENGTH] = '\0';
  //pack_kmer(kmer,packed_kmer);
  //unpack_kmer(packed_kmer,unpacked_kmer);
  //cout << kmer << "\t" << packed_kmer << "\t" << unpacked_kmer << endl;
  return 0;
}

int kad_dump(int argc, char **argv)
{
  rocksdb::DB* db;
  rocksdb::Options options;
  options.create_if_missing = true;
  rocksdb::Status status = rocksdb::DB::Open(options, LEVELDB_PATH, &db);
  assert(status.ok());

  rocksdb::Iterator* it = db->NewIterator(rocksdb::ReadOptions());
  for (it->SeekToFirst(); it->Valid(); it->Next()) {
    uint64_t kmer_int_found;
    memcpy(&kmer_int_found, it->key().data(), sizeof(uint64_t));
    std::cout << int_to_str(kmer_int_found) << ": "  << it->value().ToString() << std::endl;
  }
  return 0;
}

int kad_query(int argc, char **argv) {
  rocksdb::DB* db;
  rocksdb::Options options;
  options.create_if_missing = true;
  rocksdb::Status status = rocksdb::DB::Open(options, LEVELDB_PATH, &db);
  assert(status.ok());
  std::string value;

  if (argc < 2) {
		fprintf(stderr, "\n");
		fprintf(stderr, "Usage:   kad query [options] kmer\n\n");
		return 1;
  }

  char *kmer = argv[1];

  uint64_t kmer_int = str_to_int(kmer);

  rocksdb::Slice key((char*)&kmer_int, sizeof(uint64_t));

  rocksdb::Status s = db->Get(rocksdb::ReadOptions(), key, &value);
  if (s.ok()) {
    std::cout << kmer << "\t" << value << std::endl;
  }

  return 0;
}

int kad_index(int argc, char **argv)
{
  rocksdb::DB* db;
  rocksdb::Options options;
  options.create_if_missing = true;
  rocksdb::Status status = rocksdb::DB::Open(options, LEVELDB_PATH, &db);
  assert(status.ok());
  std::string value;

  if (argc < 2) {
		fprintf(stderr, "\n");
		fprintf(stderr, "Usage:   kad index [options] sample_name counts.tsv\n\n");
		return 1;
  }

  char *sample_name = argv[1];
  char *file = argv[2];

  gzFile fp;
	kstream_t *ks;
	kstring_t *str,*kmer;
  int dret, count;
  uint64_t kmer_int;

  kmer  = (kstring_t*)calloc(1, sizeof(kstring_t));
  str   = (kstring_t*)calloc(1, sizeof(kstring_t));
  fp = gzopen(file, "r");
  if(!fp) { fprintf(stderr, "Failed to open %s\n", file); exit(EXIT_FAILURE); }

  ks = ks_init(fp);

  while (ks_getuntil(ks, 0, str, &dret) >= 0) {
    kputs(str->s,kmer);
    if(dret != '\n') {
      if(ks_getuntil(ks, 0, str, &dret) > 0 && isdigit(str->s[0])) {
        count = atoi(str->s);
        
        kmer_int = str_to_int(kmer->s);

        rocksdb::Slice key((char*)&kmer_int, sizeof(uint64_t));

        json_t *obj; 
        json_t *counts;
        json_t *count_entry = json_object();
        json_t *count_value = json_integer(count);
        json_object_set_new(count_entry,sample_name,count_value);

        rocksdb::Status s = db->Get(rocksdb::ReadOptions(), key, &value);

        if(s.ok()) {
          obj = json_loads(value.c_str(),JSON_DECODE_ANY,NULL);
          counts = json_object_get(obj,"counts"); 
        } else {
          obj = json_object(); 
          counts = json_array();
          json_object_set_new(obj,"counts",counts);
        }

        json_array_append_new(counts,count_entry);

        char *json = json_dumps(obj,JSON_COMPACT);
        
        s = db->Put(rocksdb::WriteOptions(), key, json);

        json_decref(counts);
        json_decref(obj);
        free(json);
      }
    }
    kmer->l = 0;
  }
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
	fprintf(stderr, "         dump      Dump the KAD database\n");
	fprintf(stderr, "\n");
	return 1;
}

int main(int argc, char *argv[])
{
	if (argc == 1) return usage();
	if (strcmp(argv[1], "index") == 0) kad_index(argc-1, argv+1);
  else if (strcmp(argv[1], "dump") == 0) kad_dump(argc-1, argv+1);
  else if (strcmp(argv[1], "query") == 0) kad_query(argc-1, argv+1);
  else if (strcmp(argv[1], "test") == 0) kad_test(argc-1, argv+1);
	else {
		fprintf(stderr, "[main] unrecognized command '%s'. Abort!\n", argv[1]);
		return 1;
	}
	return 0;
}
