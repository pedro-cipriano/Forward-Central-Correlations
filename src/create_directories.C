// Pedro Cipriano, Oct 2012
// DESY, CMS
// Last Update: 13 Nov 2012
//
// create_directories(base_dir,name,type)
// create directories needed to store the results if they didn't exist yet

#include <iostream>
#include <string>
#include <sys/stat.h>
#include <sys/types.h>

using namespace std;

void create_directories(string base, string name, string mode)
{
string dir = base + name;
if (mode == "root")
{
	mkdir( dir.c_str() , S_IRWXU);
}
if (mode == "plots")
{
	mkdir( dir.c_str() , S_IRWXU);
	string dir_c = dir + "/c/";
	string dir_eps = dir + "/eps/";
	string dir_png = dir + "/png/";
	mkdir( dir_c.c_str() , S_IRWXU);
	mkdir( dir_eps.c_str() , S_IRWXU);
	mkdir( dir_png.c_str() , S_IRWXU);
}

}
