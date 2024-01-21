#ifndef UTIL_FILE_H_
#define UTIL_FILE_H_

#include <string>
#include <set>
#include <vector>
#include <iostream>
#include <fstream>
#include <cstring>

#if defined _MSC_VER
#include <direct.h>
#elif defined __GNUC__
#include <sys/stat.h>
#include <sys/types.h>
#endif

#if defined(_WIN32)
#include <direct.h>
#include <io.h>
#else
#include <dirent.h>
#include <unistd.h>
#endif

namespace util{

int GetFiles(const std::string& cate_dir,
          std::set<std::string>& files) {

#if defined(_WIN32)
  _finddata_t file;
  long lf;
  //输入文件夹路径
  if ((lf=_findfirst(cate_dir.c_str(), &file)) == -1) {
    cout<<cate_dir<<" not found!!!"<<endl;
  } else {
    while(_findnext(lf, &file) == 0) {
      //输出文件名
      //cout<<file.name<<endl;
      if (strcmp(file.name, ".") == 0 || strcmp(file.name, "..") == 0)
        continue;
      files.emplace(file.name);
    }
  }
  _findclose(lf);
#else
  DIR *dir;
  struct dirent *ptr;
  char base[1000];

  if ((dir=opendir(cate_dir.c_str())) == NULL){
    std::string error_info = std::string("Open dir error...: ") + cate_dir;
    perror(error_info.c_str());
    exit(1);
  }

  while ((ptr=readdir(dir)) != NULL){
    if(strcmp(ptr->d_name,".")==0 
    || strcmp(ptr->d_name,"..")==0)    ///current dir OR parrent dir
            continue;
    else if(ptr->d_type == 8
          || ptr->d_type == 0)    ///file
      //printf("d_name:%s/%s\n",basePath,ptr->d_name);
      files.emplace(ptr->d_name);
    else if(ptr->d_type == 10)    ///link file
      //printf("d_name:%s/%s\n",basePath,ptr->d_name);
      continue;
    else if(ptr->d_type == 4){   ///dir
      files.emplace(ptr->d_name);
      /*
            memset(base,'\0',sizeof(base));
            strcpy(base,basePath);
            strcat(base,"/");
            strcat(base,ptr->d_nSame);
            readFileList(base);
      */
    }
  }
  closedir(dir);
#endif
  return 0;
}

inline bool ExistsFile(const std::string& name) {
  std::ifstream f(name.c_str());
  if (!f.good())
    return false; 
  return true;
}

inline void Mkdir(const std::string& file_name) {
  #if defined(_WIN32)
  _mkdir(file_name.c_str());
  #else
   mkdir(file_name.c_str(),
         0777);  // notice that 777 is different than 0777
  #endif
  return;
}
  
};

#endif // UTIL_FILE_H_