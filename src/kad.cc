#include <zlib.h>
#include <iostream>
#include <rocksdb/db.h>
#include <cassert>
#include <stdlib.h>
#include <jansson.h>

#include "kseq.h"
#include "kstring.h"

KSEQ_INIT(gzFile, gzread)

#define KAD_VERSION "0.0.2"
#define LEVELDB_PATH "/tmp/testdb"

using namespace std;

int kad_dump(int argc, char **argv)
{
  rocksdb::DB* db;
  rocksdb::Options options;
  options.create_if_missing = true;
  rocksdb::Status status = rocksdb::DB::Open(options, LEVELDB_PATH, &db);
  assert(status.ok());

  rocksdb::Iterator* it = db->NewIterator(rocksdb::ReadOptions());
  for (it->SeekToFirst(); it->Valid(); it->Next()) {
    std::cout << it->key().ToString() << ": "  << it->value().ToString() << std::endl;
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

  rocksdb::Status s = db->Get(rocksdb::ReadOptions(), kmer, &value);
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
        //std::cerr << kmer->s << "\t" << count << std::endl;

        json_t *obj; 
        json_t *counts;
        json_t *count_entry = json_object();
        json_t *count_value = json_integer(count);
        json_object_set_new(count_entry,sample_name,count_value);

        rocksdb::Status s = db->Get(rocksdb::ReadOptions(), kmer->s, &value);

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
        
        s = db->Put(rocksdb::WriteOptions(), kmer->s, json);

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
	else {
		fprintf(stderr, "[main] unrecognized command '%s'. Abort!\n", argv[1]);
		return 1;
	}
	return 0;
}
