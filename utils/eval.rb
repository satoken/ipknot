#!/usr/bin/env ruby

def read_bpseq(f)
  bpseq=[]
  len=0
  IO.foreach(f) do |l|
    next if l=~/^#/
    l=l.chomp.split
    i,j=l[0].to_i,l[2].to_i
    len=i
    next if j==0
    bpseq.push([i,j]) if i<j
  end
  [bpseq,len]
end

def get_acc(tp, tn, fp, fn)
  ppv=sen=fval=mcc=0.0
  ppv=tp.to_f/(tp+fp) if tp+fp>0
  sen=tp.to_f/(tp+fn) if tp+fn>0
  fval=2*ppv*sen/(ppv+sen) if ppv+sen>0
  mcc=(tp.to_f*tn-fp*fn)/Math.sqrt((tp.to_f+fp)*(tp+fn)*(tn+fp)*(tn+fn)) if tp+fp>0 && tp+fn>0 && tn+fp>0 && tn+fn>0
  [sen, ppv, mcc, fval]
end

ans,len=read_bpseq(ARGV[0])
res,len=read_bpseq(ARGV[1])
tp=(ans&res).size
fp=res.size-tp
fn=ans.size-tp
tn=len*(len-1)/2-fp
acc=get_acc(tp,tn,fp,fn).map{|v| sprintf("%5.4f", v)}
puts (acc+[tp,tn,fp,fn,len]).join(",")
