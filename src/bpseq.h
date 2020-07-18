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

public:
  enum { 
    U = -1,   // unpaired
    DOT = -2, // . 
    L = -3,   // <
    R = -4,   // > 
    LR = -5   // |
  };

private:
  std::string seq_;
  std::vector<int> bp_;
};

#endif  //  __INC_BPSEQ_H__

