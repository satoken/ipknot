// BPSEQ
#ifndef __INC_BPSEQ_H__
#define __INC_BPSEQ_H__

#include <string>
#include <vector>

class BPSEQ
{
public:
  BPSEQ();
  BPSEQ(const std::string& seq, const std::string& paren);

  bool load(const char* fname);
  bool save(const char* fname) const;

  const std::string& seq() const { return seq_; }
  const std::vector<int>& bp() const { return bp_; }

private:
  std::string seq_;
  std::vector<int> bp_;
};

#endif  //  __INC_BPSEQ_H__

