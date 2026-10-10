#define _GNU_SOURCE
#include <sys/types.h>
#include <sys/wait.h>
#include <sys/resource.h>
#include <signal.h>
#include <unistd.h>
#include <time.h>
#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
static volatile sig_atomic_t child_pid=0, timed_out=0;
static void timeout_handler(int sig) { (void)sig; timed_out=1; if(child_pid>0) kill(-child_pid,SIGKILL); }
static double seconds(struct timespec t) { return t.tv_sec+t.tv_nsec/1e9; }
int main(int argc,char **argv) {
  if(argc<5) return 2;
  FILE *out=fopen(argv[1],"w"); if(!out) return 2;
  unsigned timeout=(unsigned)strtoul(argv[2],NULL,10);
  unsigned long memory_gib=strtoul(argv[3],NULL,10);
  struct timespec start,end; struct rusage usage; int status;
  struct sigaction sa={0}; sa.sa_handler=timeout_handler; sigemptyset(&sa.sa_mask); sigaction(SIGALRM,&sa,NULL);
  clock_gettime(CLOCK_MONOTONIC,&start);
  pid_t pid=fork(); if(pid<0) return 2;
  if(!pid) {
    setsid(); fclose(out);
    if(memory_gib) {
      struct rlimit limit={memory_gib*1024UL*1024UL*1024UL,memory_gib*1024UL*1024UL*1024UL};
      if(setrlimit(RLIMIT_AS,&limit)) { perror("setrlimit"); _exit(126); }
    }
    struct rlimit core={0,0}; setrlimit(RLIMIT_CORE,&core);
    execvp(argv[4],argv+4); perror("execvp"); _exit(127);
  }
  child_pid=pid; alarm(timeout);
  while(wait4(pid,&status,0,&usage)<0) if(errno!=EINTR) return 2;
  alarm(0); clock_gettime(CLOCK_MONOTONIC,&end);
  int code=WIFEXITED(status)?WEXITSTATUS(status):128+WTERMSIG(status);
  fprintf(out,"{\"wall_seconds\":%.9f,\"user_seconds\":%.6f,\"system_seconds\":%.6f,\"rss_kib\":%ld,\"exit_code\":%d,\"timeout\":%s}\n",seconds(end)-seconds(start),usage.ru_utime.tv_sec+usage.ru_utime.tv_usec/1e6,usage.ru_stime.tv_sec+usage.ru_stime.tv_usec/1e6,usage.ru_maxrss,code,timed_out?"true":"false");
  fclose(out); return code;
}
