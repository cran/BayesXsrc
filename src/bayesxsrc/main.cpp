#if !defined (__BUILDING_GNU)
#define __BUILDING_GNU
#endif

#if !defined (TEMPL_INCL_DEF)
#define TEMPL_INCL_DEF
#endif

#if !defined (_MSC_VER2)
#define _MSC_VER2
#endif

#if !defined (NO_TEMPLATE_FRIENDS)
#define NO_TEMPLATE_FRIENDS
#endif

#include "clstring.h"
#include "adminparse_gnu.h"
#include <iostream>
#include <fstream>
#include <string>
#include <dirent.h>
#include <sys/stat.h>
#include <unistd.h>
#include <errno.h>

#if defined(__BUILDING_LINUX)
#include <stdio.h>
#include <stdlib.h>

#if !defined (BUILD_FOR_BAYESXSRC)
#include <readline/readline.h>
#include <readline/history.h>
#endif

#endif

int main(int argc, char *argv[])
  {
  // terminating commands
//  ST::string* stop1 = new ST::string("quit") ;
//  ST::string* stop2 = new ST::string("exit") ;

  bool run=false;
  admin_gnu a;

  char *ptr;
  char path[100] = "";
  ptr = getcwd(path, 100);

  srand((unsigned)time(NULL));

#if defined(__BUILDING_LINUX)
  ST::string tempstring = ST::string(path) + "/temp";
#else
  ST::string tempstring = ST::string(path) + "\\temp";
#endif
  const char* pathtemp = tempstring.strtochar();
  int testtemp = access(pathtemp, 06);
  if(testtemp==-1)
    {
    if(errno==ENOENT)
      {
#if defined(__BUILDING_LINUX)
      mkdir(pathtemp, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
#else
      mkdir(pathtemp);
#endif
      std::cout << "NOTE: created directory " << pathtemp << endl;
      }
    else if(errno==EACCES)
      {
      std::cout << "ERROR: no write access to " << pathtemp << "!" << endl;
      return(0);
      }
    }

#if defined(__BUILDING_LINUX)
  ST::string outputstring = ST::string(path) + "/output";
#else
  ST::string outputstring = ST::string(path) + "\\output";
#endif
  const char* pathoutput = outputstring.strtochar();
  int testoutput = access(pathoutput, 00);
  if(testoutput==-1)
    {
    if(errno==ENOENT)
      {
#if defined(__BUILDING_LINUX)
      mkdir(pathoutput, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
#else
      mkdir(pathoutput);
#endif
     std::cout << "NOTE: created directory " << pathoutput << endl;
      }
    else if(errno==EACCES)
      {
      std::cout << "ERROR: no write access to " << pathoutput << "!" << endl;
      return(0);
      }
    }


//  struct stat s;
//  char pathtemp2 = char("c:\bayesx\main.cpp");
//  int testtemp2 = stat(path, &s);

//  std::cout << pathtemp << endl;
//  std::cout << testtemp << endl;

/*  DIR *pdir;
  pdir=opendir(".");
  struct dirent *pent;
  while ((pent=readdir(pdir)))
    std::cout << pent->d_name << endl;*/

//  for (int i = 0; i < argc; i++)
//    std::cout << "i=" << i << ": " <<argv[i] << " " << endl;

  bool commandline = false;
  if(argc>1)
    {
    ST::string teststring = ST::string(argv[1]);
    ST::string s = ST::string(argv[1]);
    for(int i = 2; i<argc; i++)
      s = s + " " + ST::string(argv[i]);
    int testcl = access(s.strtochar(), 04);

    if(testcl==-1)
      {
      if(errno==ENOENT)
        {
        std::cout << "NOTE: file " << s << " does not exist!" << endl;
        }
      else if(errno==EACCES)
        {
        std::cout << "Note: no read access to " << s << "!" << endl;
        return(0);
        }
      std::cout << "      proceeding in batch mode." << endl;
      }
    else
      {
      commandline = true;
      s = "usefile " + s;
      run = a.parse(s);
      }
    }

  if(!commandline)
    {
    std::cout << "BayesX - Software for Bayesian Inference in Structured Additive Regression" << endl;
    std::cout << "Version 3.0.2 (17.07.2015)" << endl;
    while(!run)
      {
      #if defined(__BUILDING_LINUX) && !defined(BUILD_FOR_BAYESXSRC)
//      rl_bind_key('\t',rl_abort);

      char *buf;
      buf = readline("BayesX> ");
      if (buf == NULL) { // EOF
        std::cout << "exiting" << std::endl;
        exit(0);
      }

      ST::string* s=new ST::string(buf);
      run = a.parse(*s);

      if (buf[0]!=0)
         add_history(buf);
     free(buf);

     #else
      std::cout << "BayesX> ";

      char array[4096];
      std::cin.getline(array, sizeof(array), '\n');
      const char* p=array;
      ST::string* s=new ST::string(p);

      run = a.parse(*s);
      #endif

      }
    } //end while run

  return(0);
  }

