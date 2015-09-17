#include "json/json.h"
#include <stdio.h>
#include <vector>

// Reads a JSON file; returns a vector of Tbars, each with a vector of partners
class Point {
    public:
    double x,y;
    };

class Partner {
  public:
    double confidence;
    Point pt;
    int z;
    };

class TBar {
  public:
    std::string status;
    double confidence;
    Point pt;
    int z;
    std::vector<Partner> partners;
    int body_id;
    };

int main(int argc, char **argv)
{
Json::Value root;   // will contains the root value after parsing.
Json::Reader reader;
FILE *fd = fopen(argv[1], "r");
const int N=1000000;
char huge[N];
fread(huge, 1, N, fd);
bool parsingSuccessful = reader.parse(huge, root );
if ( !parsingSuccessful )
{
    // report to the user the failure and their locations in the document.
    std::cout  << "Failed to parse configuration\n"
               << reader.getFormatedErrorMessages();
    return 42;
    }
std::vector<TBar> tbs;

printf("parsed file...\n");
Json::Value md = root["metadata"];
std::string des = md["description"].asString();
int         ver = md["version"].asInt();

printf("Description: '%s', version %d\n", des.c_str(), ver );
Json::Value data = root["data"];
for(unsigned int i=0; i<data.size(); i++) {
    printf("getting entry %d\n", i);
    Json::Value v = data[i];
    // read the different things that can be in 'data'.  For now only {"T-bar","partner"} is known
    Json::Value t = v["T-bar"];
    Json::Value part = v["partners"];
    if (!t.isNull()) {
        TBar tb;
	tb.status      = t["status"].asString();
	tb.confidence  = t["confidence"].asDouble();
	tb.body_id     = t["body ID"].asInt();
	Json::Value lo = t["location"];
        tb.pt.x = lo[0u].asDouble();
        tb.pt.y = lo[1u].asDouble();
        tb.z = lo[2u].asInt();
        printf("status is %s, confidence %f, loc %.1f,%lf,%d,  %d partners\n", 
	 tb.status.c_str(), tb.confidence, tb.pt.x, tb.pt.y, tb.z, part.size() );
        for(int j=0; j<part.size(); j++) {
            Partner prt;
	    prt.confidence = part[j]["confidence"].asDouble();
	    Json::Value lo = part[j]["location"];
            prt.pt.x = lo[0u].asDouble();
            prt.pt.y = lo[1u].asDouble();
            prt.z = lo[2u].asInt();
            printf("   Partner %d, confidence %f. loc=%f %f %d\n", j, prt.confidence, prt.pt.x, prt.pt.y, prt.z);
            tb.partners.push_back(prt);
	    }
        tbs.push_back(tb);
	}
    // other annotations go here
    }
// Now transfer them to the tiles
//
// First find the markers on the specified layers
// Then find all the tiles in which they occur (may be more than one by occluded edges)
// Then transform the coordinates to coordinates of that tile (use T-Bar), partner may be outside.
//   Change layers of partners to Delta from original
// Find the tile number (we already know the layer number)
// Now write the file out in json format.
//
// now print them out in json format
printf("{\"data\": [\n");
for(int i=0; i<tbs.size(); i++) {
    printf("{\"T-bar\": {\"status\": \"%s\", \"confidence\": %f, \"body ID\": %d, \"location\": [%f, %f, %d]}, \"partners\": [\n", 
     tbs[i].status.c_str(), tbs[i].confidence, tbs[i].body_id, tbs[i].pt.x, tbs[i].pt.y, tbs[i].z);
    for(int j=0; j<tbs[i].partners.size(); j++) {
	bool last = (j == tbs[i].partners.size()-1);
        printf("   {\"confidence\": %f, \"location\": [%f, %f, %d]}%c\n", 
	 tbs[i].partners[j].confidence, tbs[i].partners[j].pt.x, tbs[i].partners[j].pt.y, tbs[i].partners[j].z,
	 last ? ' ': ',' );  // comma separated if not last
	}
    bool last = (i == tbs.size()-1);
    printf("]}%c\n", last ? ' ': ',');  // comma separated if not last
    }
printf("]}\n");  // ']' ends the list of T-bars, '}' ends the 'data' section
}
