#include <leveldb/db.h>
//#include <leveldb/options.h>
#include <cassert>

using namespace std;

int main(int argc, char **argv)
{
  leveldb::DB* db;
  leveldb::Options options;
  //options.create_if_missing = true;
  leveldb::Status status = leveldb::DB::Open(options, "/tmp/testdb", &db);
  //assert(status.ok());
  return 0;
}
