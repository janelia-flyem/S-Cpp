#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <vector>
#include <algorithm>
#include <map>
#include <set>

using namespace std;
#include "json/json.h"

#define _USESTDVECTOR_ 1
#include "nr3.h"
#include "gamma.h"
#include "incgammabeta.h"
#include "fitab.h"
#include "gaussj.h"
#include "fitlin.h"
#include "svd.h"
#include "fitsvd.h"

class overlap {
  public:
    int id1, id2;
    double area;
    overlap(int i1, int i2, double a){id1 = i1; id2 = i2; area = a;}
    bool operator <(const overlap &rhs) const {return this->id1 < rhs.id1 || (this->id1 == rhs.id1 && this->id2 < rhs.id2);}
    };

class sort_overlap {
  public:
    set<overlap>::iterator it;
    // reverse sense of comparison since we want biggest first
    bool operator <(const sort_overlap &rhs) const {return (this->it)->area > (rhs.it)->area;}
    };

class where {     // This data structure sorts all connections by (from, to, Mlayer).
  public:
    int from_id;
    int to_id;
    int mLayer;              // layer in the Medulla
    mutable int sx, sy, sz;  // sums of coordinates
    mutable int n;
    where(int from, int to, int layer){from_id = from; to_id = to; mLayer = layer; n = sx = sy = sz = 0;}
    // define a less-than operator so we can make a set of these
    bool operator < (const where &rhs) const { return this->from_id < rhs.from_id || (this->from_id == rhs.from_id && this->to_id < rhs.to_id) ||
     (this->from_id == rhs.from_id && this->to_id == rhs.to_id && this->mLayer < rhs.mLayer) ;}
    };

class CellType {
  public:
    const char *name;
    int color;
    int sc_index;  // index of entry in Medulla paper
    bool new_col;  // columnar in this analysis
    CellType(const char *a, int b, bool c=false, int mp=-1) {name = a; color = b; new_col = c; sc_index = mp;}
    };

class MeanStd {  // for computing mean and standard deviation
  public:
     MeanStd(){sum = sum2 = 0.0; n=0;}
     void Reset(){sum = sum2 = 0.0; n=0;};
     void Element(double a){sum += a; sum2 += a*a; n++;}
     void Stats(double &avg, double &std){ avg = sum/n; std = sum2/n-avg*avg < 0 ? 0.0 : sqrt(sum2/n - avg*avg);}  // could be off by rounding
     double Mean(){return sum/n;}
     double Sum(){return sum;}
     double Sum2(){return sum2;}
     double Std() {double avg = sum/n; return sum2/n-avg*avg < 0 ? 0.0 : sqrt(sum2/n - avg*avg);}  // could be off by rounding
     double Var() {double avg = sum/n; return sum2/n-avg*avg < 0 ? 0.0 : sum2/n - avg*avg;}        // could be off by rounding
     int HowMany(){return n;}
  private:
     double sum, sum2;
     long int n;  // could average more the 2B elements, so need long
     };

// general purpose sorter when we want to sort the display of an array by some value
class sort_down {
  public:
    int index;    // index of entry
    int index2;   // not always used
    double val; // strength of entry
    sort_down(){index = 0; index2 = 0; val = 0.0;}
    sort_down(int i, double d){index = i; index2 = 0; val = d;}
    sort_down(int i, int j, double d){index = i; index2 = j; val = d;}
    // operator is backwards since we want biggest first
    bool operator<(const sort_down &rhs) const {return this->val > rhs.val;}
    };

class sortname {
   public:
    int index;
    double color;     // color itself is an int, but we (sometimes) add a fraction to break sort ties.
    const char *name;
    bool operator<(const sortname &rhs) const {return this->color < rhs.color || (this->color == rhs.color && strcmp(this->name,rhs.name) < 0);}
    };

class sortconn {
  public:
    const char *name;
    int count;
    double mean, std;
    bool operator<(const sortconn &rhs) const {return this->mean*this->count > rhs.mean*rhs.count;}  // highest first
    };

// general purpose sorter when we want to sort the display of an array by some value
class sorter {
  public:
    int index;    // index of entry
    double val; // strength of entry
    sorter(){index = 0; val = 0.0;}
    sorter(int i, double d){index = i; val = d;}
    // operator is backwards since we want biggest first
    bool operator<(const sorter &rhs) const {return this->val < rhs.val;}
    };

class CellStat {
  public:
    int how_many;       // How many of this cell are there
    int used_stats;     // How many were used for statistics (not marked incomplete)
    double     in_dev;  // deviation in number of inputs
    double    out_dev;  // deviation in number of outputs
    double   area_dev;  // deviation in area
    double volume_dev;  // deviation in volume
    int         color;  // set by first table, so we can see differences.
    const char *name;
    };

class Point {
  public:
    double x,y;
    Point(){x=0.0; y=0.0;}
    Point(double a, double b){x = a; y = b;}
    bool operator<(const Point &rhs) const {return this->x < rhs.x;}
    };

class Point3d {
  public:
    double x,y,z;
    Point3d(){x=0.0; y=0.0; z=0.0;}
    Point3d(double a, double b, double c){x = a; y = b; z = c;}
    bool operator<(const Point &rhs) const {return this->x < rhs.x;}
    };

class Partner {
  public:
    double confidence;
    int body_id;
    Point pt;
    int z;
    int fake_id;  // equal to body ID if named, otherwise drawn from named cells
    bool flagged;
    bool operator<(const Partner &rhs) const {return this->body_id < rhs.body_id;}
    };

class TBar {
  public:
    std::string status;
    double confidence;
    Point pt;
    int z;
    std::vector<Partner> partners;
    int body_id;
    int image_number;
    int gid;		// group id
    short roi;          // since human defined, short should be enough
    bool convergent;
    bool flagged;
    bool operator<(const TBar &rhs) const {return this->image_number < rhs.image_number;}
    };

class NeuroTr {
  public:
    const char *name;  // base name
    vector<const char *> trans;
    vector<const char *> recep;
    const char *how_tr;   // how is this known
    const char *how_re;   // how is this known
    };

class Zoid {
    public:
	Point pts[4];
    Zoid(Point p0, Point p1, Point p2, Point p3) { pts[0] = p0; pts[1] = p1; pts[2] = p2; pts[3] = p3;}
    };

vector<NeuroTr> NeuroTransmitters;

// First column is cell type, second is color, third is index in Single-column paper, if defined.
CellType CellTypes[] = {
 CellType("PPL",   87, true, 17),
 CellType("MBON",  95, true,  5),
 CellType("KC-prime",  70    ),  // Added
 CellType("KC-c",   71, true, 21),
 CellType("KC-p",   72, true,  9),
 CellType("KC-s",   73, true,  9),
 CellType("KC-a",   74, true, 11),
 CellType("Other",   55, true,  6),   // Added
 CellType("Irr",     6, true, 13),
 CellType("alpha",   24, true    ),   // columnar, but not in old
 CellType("medium",  11, true,  2),
 CellType("Mi11",  83    ),  // Added
 CellType("Mi13",   9    ),
 CellType("Mi14",  10    ),
 CellType("Mi15",  11, true, 26),
 CellType("Mi1",   12, true,  1),
 CellType("Mi2",   13    ),
 CellType("Mi4",   14, true, 19),
 CellType("Mi9",   15, true, 20),
 CellType("Tm16",  16    ),
 CellType("Tm19",  17    ),
 CellType("Tm20",  18, true, 22),
 CellType("Tm23",  19    ),
 CellType("Tm25",  20    ),
 CellType("Tm28",  79    ), // Added
 CellType("Tm29",  78    ), // Added
 CellType("Tm1",   21, true,  8),
 CellType("Tm2",   22, true, 10),
 CellType("Tm3",   23, true,  3),
 CellType("Tm4",   24, true, 12),
 CellType("Tm5Y",  48    ), // Added
 CellType("Tm5a",25    ),
 CellType("Tm5b",53, true    ), // Added, columnar but not in old
 CellType("Tm5c",54    ), // Added
 CellType("Tm6",   26, true, 16),
 CellType("Tm8",   60    ), // Added
 CellType("Tm9",   47, true, 24), // Added
 CellType("L1",    27, true,  0),
 CellType("L2",    28, true,  7),
 CellType("L3",    29, true, 18),
 CellType("L4",    30, true, 14),
 CellType("L5",    31, true,  4),
 CellType("Dm10",  58    ), // Added
 CellType("Dm11",  73    ), // Added
 CellType("Dm12",  57    ), // Added
 CellType("Dm13",  77    ), // Added
 CellType("Dm14",  59    ), // Added
 CellType("Dm15",  32    ),
 CellType("Dm16",  80    ), // Added
 CellType("Dm17",  51    ), // Added
 CellType("Dm18",  84    ), // Added
 CellType("Dm1",   76    ), // Added
 CellType("Dm2",   33, true, 25),
 CellType("Dm3",   34    ),
 CellType("Dm4",   56    ), // Added
 CellType("Dm5",   74    ), // Added
 CellType("Dm6",   81    ), // Added
 CellType("Dm:",   75    ), // Added (was 6, changed to ':' so can map 6/14/16/17 to same
 CellType("Dm7",   66    ), // Added
 CellType("Dm8",   35, true, 23),
 CellType("Dm9",   36    ),
 CellType("Pm1",   37    ),
 CellType("Pm2",   38    ),
 CellType("Pm3",   64    ), // Added
 CellType("Pm4",   65    ), // Added
 CellType("Pm",    63    ), // Added
 CellType("TmY10", 39    ),
 CellType("TmY11", 40    ),
 CellType("TmY13", 41    ),
 CellType("TmY14", 42    ),
 CellType("TmY15", 85    ), // Added BIGGEST
 CellType("TmY3",  43    ),
 CellType("TmY5",  44, true, 15),
 CellType("TmY9",  45    ),
 CellType("Y3",    46    ),
 CellType("Lawf1", 49    ),
 CellType("Lawf2", 50    ),
 CellType("CT1",   82    ),
 CellType("Dm",    68    ),  // for unclassified cells of this general type
 CellType("Tm",    69    ),  // for unclassified cells of this general type
 CellType("out",   70    ),  // for unclassified cells of this general type
 CellType("tan",   71    ),  // for unclassified cells of this general type
 CellType("ukPm",  72    ),  // for unclassified cells of this general type
 CellType("Unknown",67   ) };

int NumTypes = sizeof(CellTypes)/sizeof(CellType);

// Translate a sequential number to a color
int clr(int n) {
int tbl4[] = {0x40, 0x80, 0xc0, 0xFF};
int tbl5[] = {0x30, 0x60, 0x90, 0xC0, 0xF0};
const int NRED   = 5;
const int NGREEN = 5;
const int NBLUE  = 4;
int red =    n % NRED;
n = n / NRED;
int green =  n % NGREEN;
n = n / NGREEN;
int blue  =  n % NBLUE;
int rslt = tbl5[red];
rslt   += (tbl5[green] << 8);
rslt   += (tbl4[blue] << 16);
return rslt;
}

int color(const char *name) {
for(int k=0; k<NumTypes; k++) {
    if (strncmp(CellTypes[k].name, name, strlen(CellTypes[k].name)) == 0) 
	return clr(CellTypes[k].color);
    }
return 0xFFFFFF;
}

// If we don't find it otherwise, it's Unknown
const char *BaseName(const char *name) {
for(int k=0; k<NumTypes-1; k++) {
    if (strncmp(CellTypes[k].name, name, strlen(CellTypes[k].name)) == 0)  {
        // if (CellTypes[k].color == 63) {  // Just testing
	    // printf("Bingo! %s\n", name);
	    // exit(42);
	    // }
	return CellTypes[k].name;
	}
    }
return "Unknown";
}

// Anything we do not find must return the last entry, Unknown.
int TypeIndex(const char *name) {
for(int k=0; k<NumTypes-1; k++) {
    if (strncmp(CellTypes[k].name, name, strlen(CellTypes[k].name)) == 0) 
	return k;
    }
return NumTypes-1;
}

const int NYELLOW = 4;  // a, c, e, f
const int NPALE   = 3;  // H, b, d

// Is it a yellow column?  It is if it ends in a, c, e, or f.
// Also need to check if the char before is not a lower case letter, otherwise 'Dm8outside' gets column e
bool IsYell(const char *name) {
int n = strlen(name) -1;
//bool OKsoFar = n-1 > 0 && (name[n-1] < 'a' || name[n-1] > 'z');
bool OKsoFar = n-1 > 0 && (name[n-1] == '-');
return OKsoFar && (name[n] == 'A' || name[n] == 'C' || name[n] == 'E' || name[n] == 'F');
}
// Is it a Home column cell?  It is if it ends in " home"
bool IsHome(const char *name) {
int n = strlen(name) - 5;
return (n > 0 && strcmp(name+n, " home") == 0);
}
// Is it a pale column?  It is if it ends in b or d, or it's a home column
// Also need to check if the char before is not a lower case letter, otherwise 'Dm8outboard' gets column d
bool IsPale(const char *name) {
int n = strlen(name) -1;
bool OKsoFar = n-1 > 0 && (name[n-1] == '-');
return ( OKsoFar && (name[n] == 'B' || name[n] == 'D')) || IsHome(name);
}
// Is it one of the 7 main columns
bool Is7col(const char *name) {
return IsPale(name) || IsYell(name);
}

// Return TRUE if its a type associated with a single column.  That means in the original table, but not T4 or Tm3
bool IsColumnar(const char *name) {
if (strncmp(name, "T4",2) == 0) return false;
if (strncmp(name, "Tm3",3) == 0) return false;
if (strncmp(name, "Tm4",3) == 0) return false;
int n = TypeIndex(name);
return CellTypes[n].new_col;
}

// Return the column associated with a cell, if this makes sense and the column is known.
// Returns a-f for surrounding cells, H for home column cell, ' ' otherwise.
char Column(const char *name) {
//if (strcmp(name, "T3") == 0)  return ' ';    // These four types, though in the matrix, are not columnnar
//if (strcmp(name, "T4") == 0)  return ' ';
if (strncmp(name, "Tm3",3) == 0 || strncmp(name, "Tm4",3) == 0) {
    char c = ' ';
    if (strlen(name) >= 5) {
        if (strncmp(name+4, "home", 4) == 0)
	    c = 'H';
        else if (name[4] >= 'A' && name[4] <= 'F')
	    c = tolower(name[4]);
	}
    printf("Tm Testing '%s' -> '%c'\n", name, c);
    return c;
    }
if (IsHome(name))             return 'H';
int n = TypeIndex(name);
bool is_c = CellTypes[n].new_col;
char c = name[strlen(name)-1];
if (is_c  && (c >= 'A' && c <= 'F')) return tolower(c);
return ' ';
}

// Returns a numeric version of the Column ID, for creating matrices.
// 0 = Unknown or not columnar, 1 = Home, 2-7 = a-f

int NumColumn(const char *name) {
char c = Column(name);
if (c == ' ') return 0;
if (c == 'H') return 1;
if (c >= 'a' && c <= 'f') return c - 'a' + 2;
printf("Bad code %c in NumColumn\n", c);
exit(42);
}

// Print stats about an entry in HTML.  If none, print 3 blanks.  If one, print count and mean.
// If more than one, add std deviation.
// If the 
void PrintEntries(FILE *fp, MeanStd m, double *all = NULL, double *weak = NULL){
char c = ' ';
if (all != NULL) { // if we want to add to statistics
    *all += m.Sum();
    if (m.HowMany() <= 3) {
	*weak += m.Sum();
        c = '+';
	}
    }
if (m.HowMany() == 0)
    fprintf(fp, "<TD><TD><TD>\n");
else if (m.HowMany() == 1)
    fprintf(fp, "<TD>%d<TD>%.2f%c<TD>\n", m.HowMany(), m.Mean(), c);
else {
    fprintf(fp, "<TD>%d<TD>%.2f%c<TD>%.2f\n", m.HowMany(), m.Mean(), c, m.Std());
    }
}

// Take a vector of doubles, return a vector of bools that indicates the 7 (at most) largest.
// Do not include any that are marked incomplete
vector<bool> Largest(vector<double> &data, vector<bool> &incomplete) {
int N = data.size();
vector<bool> rslt(N, false);
vector<sort_down> sa;
for(int i=0; i<data.size(); i++) {
    if (!incomplete[i])
	sa.push_back(sort_down(i, data[i]));
    }
sort(sa.begin(), sa.end() );
for(int i=0; i<7 && i<sa.size(); i++)
    rslt[sa[i].index] = true;
return rslt;
}

// writes into the array list all the known inputs or outputs of the neuron type.
// If we are looking at INPUTs, then we want the transmitter types of the cells.
// If we are looking at outputs, we want the receptor types.
void GetNlist(char *list, char *title, const char *name, bool input) {
strcpy(list,  ""); // assume it's blank
strcpy(title, ""); // assume it's blank
for(int i=0; i<NeuroTransmitters.size(); i++) {
    if (strcmp(name, NeuroTransmitters[i].name) != 0)
	continue;
    vector<const char *>l = input ? NeuroTransmitters[i].trans : NeuroTransmitters[i].recep;
    for(int j=0; j<l.size(); j++) {
        const char *p = strchr(l[j], '_');
        if (p == NULL)
	    strcat(list, l[j]);  // if no _, just copy
        else {
	    char local[128];
	    strcpy(local, l[j]);
	    char *p1 = strtok(local, "_");
            char *p2 = strtok(NULL,  "\n");
            strcat(list, p1); strcat(list, "<sub>"); strcat(list, p2); strcat(list,"</sub>");
	    }
        if (j < l.size()-1)
	    strcat(list, ",");
	}
    strcpy(title, input ? NeuroTransmitters[i].how_tr : NeuroTransmitters[i].how_re);
    return;
    }
}

set<where> ss;  // ss = set of synapses

double GetZ(int from, int to) {
int zsum = 0;
int N = 0;
for(int i=1; i<=10; i++) {
    where local(from, to, i);
    set<where>::iterator it = ss.find(local);
    if (it != ss.end()) {
	zsum += it->sz;
        N += it->n;
	}
    }
if (N == 0) {
    printf("Nothing in table with from=%d and to=%d\n", from, to);
    exit(42);
    }
return double(zsum)/N/100.0;  // return value in microns
}

// Convert a Z coordinate to an (approximate) layer
int ZtoLayer(double z) {
int m = 0;
if (z > 14.0) m =  1;
if (z > 20.5) m =  2;
if (z > 27.5) m =  3;
if (z > 34.5) m =  4;
if (z > 38.5) m =  5;
if (z > 42.5) m =  6;
if (z > 47.5) m =  7;
if (z > 52.5) m =  8;
if (z > 59.5) m =  9;
if (z > 69.5) m = 10;
return m;
}

vector<vector<double> > planes(11, vector<double>(3,0));

// More sophisticated layer calculation, using the planes defined by Shinya's bookmarks.
int XYZtoLayer(double x, double y, int z) {
int m = 1;
for(int k=1; k<=9; k++) {
    // planes[k] defines the Mk/Mk+1 boundary
    double zpl = planes[k][0] + planes[k][1]*x + planes[k][2]*y;
    if (z > zpl)
	m = k+1;
    }
return m;
}
// More sophisticated layer calculation, using the planes defined by Shinya's bookmarks.
// This version also allows layer 0 (below M1) and 11 (> M10)
int XYZtoLayer011(double x, double y, int z) {
int m = 0;
for(int k=0; k<=10; k++) {
    // planes[k] defines the Mk/Mk+1 boundary
    double zpl = planes[k][0] + planes[k][1]*x + planes[k][2]*y;
    if (z > zpl)
	m = k+1;
    }
return m;
}

// Returns a copy of a string with ' ' replaced by '-'.  Need to make a copy if you call it more than once.
char nbuf[128];
char *NoSpace(const char *str)
{
char *q = nbuf;
for(const char *p = str; *p != 0; p++)
    *q++ = (*p == ' ') ? '-' : *p;
*q = 0;
return nbuf;
}

vector<TBar> tbs;
// Given two body IDs, find the ROI of (the first) synapse.
// Super slow, needs to be replaced....
int GetRoi(int from, int to, map<int,int> &rev_map)
{
from = rev_map[from];  // first, convert from table indices to body IDs
to   = rev_map[to];
for(int i=0; i<tbs.size(); i++) {
    if (tbs[i].body_id != from)
	continue;
    for(int j=0; j<tbs[i].partners.size(); j++) {
	if (tbs[i].partners[j].body_id == to)
	    return tbs[i].roi;
	}
    }
return 0;
}

vector<CellStat> CellStats;

// write the html file.  Name it after the shortest of the row names.
void write_html(vector<char *> & names, map<int,int> &rev_map, vector<int>&row_ids, const char *suffix, 
 vector<vector<int> > &rows, vector< vector<int> > &dat, vector<vector<int> > &recip, 
 vector<int> &sums, vector<double> &volumes, vector<double> &areas, vector<bool> incomplete,
 vector<MeanStd> &stats, vector<MeanStd> &seven, vector<MeanStd> &pale, vector<MeanStd> &yellow,
 vector<vector<int> > &col_matrix,  /* 8x8 matrix of cross-column connections*/ vector<int> &Tbar_count)
{
if (row_ids.size() <= 0)
    return;
FILE *fp;
char fname[128];
const char *bname = BaseName(names[row_ids[0]]);
sprintf(fname,"%s-%s.html", suffix, bname);
for(char *p = fname; *p != '\0'; p++) {
    if (*p == '/') *p = '-';             // replace '/' by -
    if (*p == ' ') *p = '_';             // replace ' ' by _
    }
fp = fopen(fname,"w");
if (fp == NULL) {
    printf("could not open %s for write\n", fname);
    exit(42);
    }


int width = 0; // number of entries in longest row
for(int i=0; i<dat.size(); i++) {
    width = max(width, int(dat[i].size()));
    }
const int HOW_MANY_TO_PRINT = 100;
int print_width = min(width, HOW_MANY_TO_PRINT);  // print at most 12

int total = 0; // total strength of all entries in table
int super = 0;
for(int i=0; i<dat.size(); i++) {
    for(int k=0; k<dat[i].size() && k < print_width; k++)
        total += dat[i][k];
    super += sums[i];
    }

// This vector will tell the largest (up to 7 of them) that are complete.  Find the mean volume of the completed cells.
vector<bool> mask = Largest(volumes, incomplete);
MeanStd CompVol;
for(int i=0; i<volumes.size(); i++) {
    if (!incomplete[i])
	CompVol.Element(volumes[i]);
    }

bool input = strcmp(suffix, "in") == 0;
fprintf(fp, "<!DOCTYPE html PUBLIC \"-//W3C//DTD HTML 4.01 Transitional//EN\">\n");
fprintf(fp, "<html>\n");
fprintf(fp, "  <head>\n");
fprintf(fp, "    <meta http-equiv=\"content-type\" content=\"text/html;\n");
fprintf(fp, "      charset=windows-1252\">\n");
fprintf(fp, "    <title> %s Neurons of Type %s</title>\n", input? "Inputs to" : "Outputs from", bname);
fprintf(fp, "   <style type=\"text/css\">a {text-decoration: none; color: #000000; font-weight:bold;}</style>\n"); // non-underlined links
fprintf(fp, "  </head>\n");
fprintf(fp, "  <body>\n");
// Contents go here
set<int> colors;   // will be all colors used in the table
fprintf(fp,"<p>\n");
fprintf(fp, "<h1>%s Neurons of Type %s</h1>\n", input? "Inputs to" : "Outputs from", bname);
fprintf(fp,"</p><a href=\"%s-%s.html\"> Click here to see %s of %s instead.</a><p>\n",
 input? "out":"in", bname, input?"outputs":"inputs", bname);
fprintf(fp,"<p>Table of all cells of type %s, one per row, %d rows.  Cells with the cell name in green are thought at least 98%% complete by the annotator.  Cells in yellow have at least 90%% of the volume of those marked complete.</p>", bname, row_ids.size());
if (width <= HOW_MANY_TO_PRINT)
    fprintf(fp,"  All known %s (max %d) are shown in the table.\n", input ? "inputs" : "outputs", width);
else {
    fprintf(fp,"  The entries in the row are the top %d %s the cell, ordered by strength.\n", 
     print_width, input ? "inputs to" : "outputs from");
    fprintf(fp,"The number %d is arbitrary, to keep the table a manageable width.\n", print_width);
    fprintf(fp,"These top %d entries cover %.2f%% of all connections to cells of this type.\n", print_width, double(total)/super*100.0);
    }
fprintf(fp,"</p><p>\n");
fprintf(fp, "<TABLE BORDER=4 CELLSPACING=4 CELLPADDING=4>\n");
if (strcmp(suffix,"in") == 0)
    fprintf(fp, "<caption> Number of connections (PSD in %s, TBar in named call) : directly reciprocal count </caption>\n", bname);
else
    fprintf(fp, "<caption> Number of connections (TBar in %s, PSD in named call) : directly reciprocal count </caption>\n", bname);
fprintf(fp, "<TD> Cell type ");
for(int k=0; k<print_width; k++)
    fprintf(fp, "<TD> #%d ", k+1);
fprintf(fp, "<TD>Conn shown/total<TD>T-bars<TD>Volume<TD> Area</TD></TR>\n");
for(int i=0; i<row_ids.size(); i++) {

    // we'll use a light yellow for almost complete cells, a light green for complete, white for incomplete
    int color_row = incomplete[i] ? (CompVol.HowMany() > 0 && volumes[i] > 0.9 * CompVol.Mean() ? 0xFFFF99 : 0xFFFFFF) : 0x99FF99;
    //int color_row = incomplete[i] ? 0xA0A0A0 : 0xFFFFFF;  // grey if incomplete
    fprintf(fp, "<TR>\n");
    if (incomplete[i])
        fprintf(fp, "<TD BGCOLOR=\"#%06x\"> %s %d\n", color_row, names[row_ids[i]], rev_map[row_ids[i]] );
    else  // complete, so print in bold
        fprintf(fp, "<TD BGCOLOR=\"#%06x\"> <strong>%s</strong> %d\n", color_row, names[row_ids[i]], rev_map[row_ids[i]] );
    int shown = 0;
    for(int j=0; j<rows[i].size() && j < print_width; j++) {
        char linkname[128];
        sprintf(linkname, "%s-%s.html", suffix, BaseName(names[rows[i][j]]) );
        int from = input ? rows[i][j] : row_ids[i];
        int   to = input ? row_ids[i] : rows[i][j];
	fprintf(fp, "<TD BGCOLOR=\"#%06x\"><a href=\"%s\"><pre>%s(&alpha;%d)<br>%4d:%2d<br>%d</pre></a>\n", 
         color(names[rows[i][j]]), linkname, names[rows[i][j]], GetRoi(from, to, rev_map), dat[i][j], recip[i][j], rev_map[rows[i][j]]);
        colors.insert(color(names[rows[i][j]]));
        shown += dat[i][j];
	}
    for(int j = rows[i].size(); j<print_width; j++)
	fprintf(fp,"<TD>\n");
    fprintf(fp, "<TD BGCOLOR=\"#%06x\"><pre> %3d/%3d </pre><TD>%d<TD BGCOLOR=\"#%06x\">%.2f %c%c<TD BGCOLOR=\"#%06x\">%.2f\n", 
     color_row, shown, sums[i], Tbar_count[row_ids[i]], color_row, volumes[i], incomplete[i]?'?':' ', mask[i]?'*':' ', color_row, areas[i]);
    fprintf(fp, "</TD> </TR>\n");
    }
// Now write a final row, with the averages.
MeanStd su, vu, ar;
for(int i=0; i<row_ids.size(); i++) {
    if (mask[i]) {
        su.Element(sums[i]);
        vu.Element(volumes[i]);
        ar.Element(areas[i]);
        }
    }
fprintf(fp, "<TD> Cell type ");
for(int k=0; k<print_width; k++)
    fprintf(fp, "<TD> #%d ", k+1);
fprintf(fp, "<TD>%.1f&plusmn%.1f<br>%.1f%%<TD>%.1f&plusmn;%.1f<br>%.1f%%<TD>%.1f&plusmn;%.1f<br>%.1f%%</TD></TR>\n",
 su.Mean(), su.Std(), su.Std()/su.Mean()*100, vu.Mean(), vu.Std(), vu.Std()/vu.Mean()*100, ar.Mean(), ar.Std(), ar.Std()/ar.Mean()*100 );
fprintf(fp, "</TABLE>\n");
fprintf(fp,"</p>\n");

//Now re-sort each row by Z.  First expand each entry into a set of iterators to layer specific connections
vector<vector<set<where>::iterator> > sz;
for(int i=0; i<row_ids.size(); i++) {
    vector<set<where>::iterator> s;
    for(int j=0; j<rows[i].size(); j++) {
        int from = input ? rows[i][j] : row_ids[i];
        int   to = input ? row_ids[i] : rows[i][j];
        for(int k=1; k<=10; k++) {
            where local(from, to, k);
            set<where>::iterator it = ss.find(local);
            if (it != ss.end())
		s.push_back(it);
	    }
	}
    printf("Expansion is %d to %d\n", rows[i].size(), s.size() );
    sz.push_back(s);
    }
// At this point, all the entries s[i][*] will have the same from (if doing outputs) or the same 'to'.

// first delete any entry with strength <= THRESHOLD (they clutter the table too much, and hard to indentify
// if that is still not enough, try with a bigger threshold
const int THRESHOLD=2;
int thresh = THRESHOLD;
int mx = 0;
do {
    for(int i=0; i<row_ids.size(); i++) {
	for(int j = sz[i].size()-1; j >= 0; j--) {
            where look(sz[i][j]->to_id, sz[i][j]->from_id, sz[i][j]->mLayer);
            set<where>::iterator it = ss.find(look);
            int recip = it != ss.end() ? it->n : 0;
	    if (sz[i][j]->n <= thresh && recip <= thresh) { // zap it
		sz[i].erase(sz[i].begin()+j);
		}
	    }
	}
    mx = 0;
    for(int i=0; i<sz.size(); i++)
	mx = max(mx, int(sz[i].size()));
    thresh++;
    }
while (mx > 200);

// Now sort each row by the Z of its entries
vector<vector<sorter> > save;
for(int i=0; i<row_ids.size(); i++) {
    vector<sorter> sa(sz[i].size());
    for(int j=0; j<sz[i].size(); j++) {
        sa[j].val = double(sz[i][j]->sz)/sz[i][j]->n/100.0;  // convert to microns
        sa[j].index = j;
	}
    sort(sa.begin(), sa.end());
    save.push_back(sa);
    vector<set<where>::iterator> copy_sz = sz[i];
    for(int j=0; j<sz[i].size(); j++) {
	sz[i][j] = copy_sz [sa[j].index];
        }
    }
// approximately align each table
double margin = strcmp(bname, "Mi1") == 0 ? 1.6 : 1.0;  // use a bigger value for Mi1 to avoid end of table
for(int i=0; i<1000; i++) {
    // for column i, find the number of elements and the smallest z
    int num = 0;
    double smallz = 1e30;
    for(int j=0; j<row_ids.size(); j++) {
	if (sz[j].size() > i) { // data exists
	    num++;
            smallz = min(smallz, save[j][i].val);
	    }
	}
    if (num < 2) break;
    printf("Column %d, %d are used, smallest %.2f\n", i, num, smallz);
    // now, for every one that exists, and is more than smallz+margin, bump it over.
    for(int j=0; j<row_ids.size(); j++) {
	if (sz[j].size() > i) { // data exists
            if(save[j][i].val > smallz+margin) {
                sz  [j].insert(  sz[j].begin()+i, sz[j][0]     );  // Need to insert SOME element; will never be referenced
                save[j].insert(save[j].begin()+i, sorter(0,0.0));  // since we will check vector 'save' before de-referencing
		}
	    }
	}
    }
// re-compute how many to print
width = 0; // number of entries in longest row
for(int i=0; i<dat.size(); i++) {
    width = max(width, int(sz[i].size()));
    }
print_width = min(width, HOW_MANY_TO_PRINT);  // print at most HOW_MANY_TO_PRINT

// Find an average Z for each column
vector<double> Zavg(print_width);
for(int j=0; j<print_width; j++) {
    MeanStd m;
    for(int i=0; i<row_ids.size(); i++) {
	if (sz[i].size() > j && save[i][j].val > 0.0) {
	    m.Element(double(sz[i][j]->sz)/sz[i][j]->n/100.0);  // convert to microns
	    }
	}
    Zavg[j] = m.Mean();
    }
// Find an approximate Medulla layer
vector<int> Mlayer(print_width);
for(int i=0; i<print_width; i++) {
    Mlayer[i] = ZtoLayer(Zavg[i]);;
    }

// Reprint the table, now sorted by Z
fprintf(fp,"<p>\n");
fprintf(fp,"Table of all cells of type %s, one per row.  The entries in the row are the top %d %s the cell, ordered by Z.\n", 
 bname, print_width, input ? "inputs to" : "outputs from");
fprintf(fp,"The number %d is arbitrary, to keep the table a manageable width.\n", print_width);
fprintf(fp,"These top %d entries cover %.2f%% of all connections to cells of this type.\n", print_width, double(total)/super*100.0);
fprintf(fp,"</p><p>\n");
fprintf(fp, "<TABLE BORDER=4 CELLSPACING=4 CELLPADDING=4>\n");
if (strcmp(suffix,"in") == 0)
    fprintf(fp, "<caption> Number of connections (PSD in %s, TBar in named call) : directly reciprocal count </caption>\n", bname);
else
    fprintf(fp, "<caption> Number of connections (TBar in %s, PSD in named call) : directly reciprocal count </caption>\n", bname);
fprintf(fp, "<TD> Cell type ");
for(int k=0; k<print_width; k++) {
    int color = (Mlayer[k] & 1) ? 0xFFFFFF : 0xCCCCCC;
    fprintf(fp, "<TD BGCOLOR=\"#%06x\"> z=%.1f <br> M%d", color, Zavg[k], Mlayer[k]);
    }
fprintf(fp, "<TD>Conn shown/total<TD>Volume<TD> Area</TD></TR>\n");
for(int i=0; i<row_ids.size(); i++) {

    // we'll use a light yellow for almost complete cells, a light green for complete, white for incomplete
    int color_row = incomplete[i] ? (CompVol.HowMany() > 0 && volumes[i] > 0.9 * CompVol.Mean() ? 0xFFFF99 : 0xFFFFFF) : 0x99FF99;
    fprintf(fp, "<TR>\n");
    fprintf(fp, "<TD BGCOLOR=\"#%06x\"> %s\n", color_row, names[row_ids[i]]);
    int shown = 0;
    for(int j=0; j<sz[i].size() && j < print_width; j++) {
        if (save[i][j].val <= 0.0) {
	    fprintf(fp,"<TD>\n");
	    continue;
            }
        int entry = input ? sz[i][j]->from_id : sz[i][j]->to_id;
        char linkname[128];
        sprintf(linkname, "%s-%s.html", suffix, BaseName(names[entry]) );
        // Find the strength of the reciprocal
	where look(sz[i][j]->to_id, sz[i][j]->from_id, sz[i][j]->mLayer);
	set<where>::iterator it = ss.find(look);
	int recip = it != ss.end() ? it->n : 0;
	//fprintf(fp, "<TD BGCOLOR=\"#%06x\"><a href=\"%s\"><pre>%s(%.1f)<br>%4d:%2d</pre></a>\n", 
         //color(names[rows[i][j]]), linkname, names[rows[i][j]], GetZ(from, to), dat[i][j], recip[i][j]);
	fprintf(fp, "<TD BGCOLOR=\"#%06x\"><a href=\"%s\"><pre>%s<br>%4d:%2d</pre></a>\n", 
         color(names[entry]), linkname, names[entry], sz[i][j]->n, recip);
        colors.insert(color(names[entry]));
        shown += sz[i][j]->n;
	}
    for(int j = sz[i].size(); j<print_width; j++)
	fprintf(fp,"<TD>\n");
    // we'll use a light yellow for incomplete, a light green for complete
    fprintf(fp, "<TD BGCOLOR=\"#%06x\"><pre> %3d/%3d </pre><TD BGCOLOR=\"#%06x\">%.2f %c%c<TD BGCOLOR=\"#%06x\">%.2f\n", 
     color_row, shown, sums[i], color_row, volumes[i], incomplete[i]?'?':' ', mask[i]?'*':' ', color_row, areas[i]);
    fprintf(fp, "</TD> </TR>\n");
    }
// Now write a final row, with the averages.
fprintf(fp, "<TD> Cell type ");
for(int k=0; k<print_width; k++)
    fprintf(fp, "<TD> #%d ", k+1);
fprintf(fp, "<TD>%.1f&plusmn%.1f<br>%.1f%%<TD>%.1f&plusmn;%.1f<br>%.1f%%<TD>%.1f&plusmn;%.1f<br>%.1f%%</TD></TR>\n",
 su.Mean(), su.Std(), su.Std()/su.Mean()*100, vu.Mean(), vu.Std(), vu.Std()/vu.Mean()*100, ar.Mean(), ar.Std(), ar.Std()/ar.Mean()*100 );
fprintf(fp, "</TABLE>\n");
fprintf(fp,"</p>\n");
    
// Add to vector used to sort by deviation(s).  Input gets called first, so it allocates.  Output adds to the end record.
if(input) {
    CellStat s;
    s.how_many = row_ids.size();
    s.used_stats = su.HowMany();
    s.name = bname;
    // Since we drop incomplete cells, it's possible there are none left.  Set deviations to 0 in this case, so we still can sort
    // without the complications of NaNs.
    if (su.HowMany() > 0) {
        s.in_dev =     su.Std()/su.Mean();
        s.volume_dev = vu.Std()/vu.Mean();
        s.area_dev   = ar.Std()/ar.Mean();
	}
    else
	s.in_dev = s.volume_dev = s.area_dev = 0.0;
    CellStats.push_back(s);
    }
else {
    if (su.HowMany() > 0)
        CellStats[CellStats.size()-1].out_dev = su.Std()/su.Mean();
    else
        CellStats[CellStats.size()-1].out_dev = 0.0;
    }

// Now for a table of cross-column connections.
int summary[] = {
// U  H  a  b  c  d  e  f
   0, 0, 0, 0, 0, 0, 0, 0,  // U
   0, 0, 1, 2, 3, 4, 5, 6,  // H
   0, 4, 0, 3, 0, 0, 0, 5,  // a
   0, 5, 6, 0, 4, 0, 0, 0,  // b
   0, 6, 0, 1, 0, 5, 0, 0,  // c
   0, 1, 0, 0, 2, 0, 6, 0,  // d
   0, 2, 0, 0, 0, 3, 0, 1,  // e
   0, 3, 2, 0, 0, 0, 4, 0   // f
   };
int scolors[] = {0xFFFFFF, 0xFFFF88, 0xFF88FF, 0x88FFFF, 0xFF8888, 0x88FF88, 0x8888FF};

vector<int> dirs(7,0);

const char *labels = "UHabcdef";
fprintf(fp,"<p><pre>\n");
fprintf(fp,"       Posterior\n");
fprintf(fp,"\n");
fprintf(fp,"          b c\n");
fprintf(fp,"Dorsal   a H d   Ventral\n");
fprintf(fp,"          f e\n");
fprintf(fp,"\n");
fprintf(fp,"       Anterior\n");
fprintf(fp,"</pre>\n");

fprintf(fp,"Table of cross column connections.  U = unknown or not columnar, H = home column, a-f surrounding columns.  Connections between nearest neighbors are indicated in color - one color for each of the 6 directions possible for a nearest neighbor.\n");
fprintf(fp,"<TABLE> <TR> <TD>\n"); // giant table containing two tables.
fprintf(fp, "<TABLE BORDER=4 CELLSPACING=4 CELLPADDING=4>\n");
fprintf(fp, "<colgroup>\n");
fprintf(fp, "    <col style=\"width: 20%\" />\n");
for(int k=0; k<8; k++)
    fprintf(fp, "    <col style=\"width: 10%\" />\n");
fprintf(fp, "  </colgroup>\n");

fprintf(fp, "<caption> Cross column connections </caption>\n");
fprintf(fp, "<TR><TD> Column ");
for(int k=0; k<8; k++)
    fprintf(fp, "<TD> %c ", labels[k]);
fprintf(fp, "</TD></TR>\n");
for(int i=0; i < 8; i++) {
    fprintf(fp, "<TR>\n");
    fprintf(fp, "<TD> %s %c", input ? "To" : "From", labels[i]);
    for(int j=0; j<8; j++) {
        int dir = summary[i*8+j];
        dirs[dir] += col_matrix[i][j];
        int color = (i == j) ? 0xE0E0E0 : scolors[dir];
        if (col_matrix[i][j] == 0)
	    fprintf(fp, "<TD BGCOLOR=\"#%06x\"> ", color );
        else
	    fprintf(fp, "<TD BGCOLOR=\"#%06x\">%d", color, col_matrix[i][j]);
	}
    fprintf(fp, "</TD> </TR>\n");
    }
fprintf(fp, "</TABLE>\n");
fprintf(fp, "<TD> <p>Each of the 6 basic directions occurs 4 times in our 7-column data.  For example, a->H, H->d, b->c, and f->e all represent the same displacement vector.  Here we sum the 4 entries that have the same direction, in each of the 6 basic directions.</p>\n");

// Table of summary by direction.
fprintf(fp, "<TABLE BORDER=4 CELLSPACING=4 CELLPADDING=4>\n");
fprintf(fp, "<caption> Summary of connections by direction</caption>\n");
for(int k=1; k<=6; k++) {
    if (input)
        fprintf(fp, "<TR><TD BGCOLOR=\"#%06x\"> %c->H <TD> %d </TD></TR>\n", scolors[k], 'a'+k-1, dirs[k]);
    else
        fprintf(fp, "<TR><TD BGCOLOR=\"#%06x\"> H->%c <TD> %d </TD></TR>\n", scolors[k], 'a'+k-1, dirs[k]);
    }
fprintf(fp, "</TABLE>\n");

fprintf(fp, "</TD></TR>\n");  // end the giant table
fprintf(fp, "</TABLE>\n");

fprintf(fp,"</p>\n");

// Now write another table, with each color, a cell of that color, and some data about it
// This table (a sorted version of info from the first table) seems not so useful.  So ignore it for now.
if (false) {
    vector<sortconn>s;
    set<int>::iterator it;
    for(it = colors.begin(); it != colors.end(); it++) {
	// go through the table, collecting data about that color
	MeanStd m;
	const char *example;
	for(int i=0; i<rows.size(); i++) {
	    for(int j=0; j<rows[i].size(); j++) {
		if (color(names[rows[i][j]]) == *it) {
		    example = names[rows[i][j]];
		    m.Element(dat[i][j]);
		    }
		}
	    }
	sortconn ss;
	ss.name = example;
	ss.count = m.HowMany();
	m.Stats(ss.mean, ss.std);
	s.push_back(ss);
	}
    sort(s.begin(), s.end());   // sort by mean strength

    fprintf(fp,"<p>\n");
    fprintf(fp, "<TABLE BORDER=4 CELLSPACING=4 CELLPADDING=4>\n");
    fprintf(fp, "<caption> Statistics of the connections between cell types, from table </caption>\n");
    fprintf(fp, "<TD> Cell type <TD> Count <TD> Per inst. <TD> Mean <TD> Std. Dev <TD> Dev as %% of Mean<TD> Total %%"
    "</TD></TR>\n");
    for(int i=0; i<s.size(); i++) {
        const char *p = BaseName(s[i].name);
        fprintf(fp, "<TD BGCOLOR=\"#%06x\"><pre>%s</pre>\n", color(s[i].name), p);
        fprintf(fp, "<TD> %d <TD> %.2f <TD> %.2f\n", s[i].count, double(s[i].count)/row_ids.size(), s[i].mean); //mean
        if (s[i].count == 1)
	    fprintf(fp,"<TD> <TD>");  // no std if only one.
        else {
            fprintf(fp, "<TD> %.2f\n", s[i].std); //std
            fprintf(fp, "<TD> %.2f%%\n", s[i].std/s[i].mean*100.0);  //frac
            }
        fprintf(fp, "<TD> %.2f%%\n", s[i].mean*s[i].count/double(total)*100.0);
        fprintf(fp, "</TD> </TR>\n");
        };
    fprintf(fp, "</TABLE>\n");
    fprintf(fp,"</p>\n");
    }
//-------------------------------------- Now write the table of ALL connections, not just those in table
total = 0;                          // Now re-do total to include all
for(int i=0; i<NumTypes; i++) {
    if (stats[i].HowMany() > 0) {
        double mean, std;
        stats[i].Stats(mean, std);
        total += int(stats[i].HowMany() * mean + 0.5);  // should be int, but round...
        printf(" Stats %d, %f\n", stats[i].HowMany(), mean);
	}
    }
vector<sort_down> sa(NumTypes);
for(int i=0; i<NumTypes; i++) {
    sa[i].index = i;
    if (stats[i].HowMany() > 0)
        sa[i].val = stats[i].Mean()*double(stats[i].HowMany())/row_ids.size();
    else
        sa[i].val = 1.0;  // don't care how these sort, since they are not printed, but must be well defined.
    }
for(int i=0; i<NumTypes; i++)
    printf("sa[%d] %s %s %d %f\n", i, bname, suffix, sa[i].index, sa[i].val);
sort(sa.begin(), sa.end()-1);
for(int i=0; i<NumTypes; i++)
    printf("sa[%d] %s %s %d %f\n", i, bname, suffix, sa[i].index, sa[i].val);
fprintf(fp,"<p>\n");
fprintf(fp, "<TABLE BORDER=4 CELLSPACING=4 CELLPADDING=4>\n");
char list[128], title[256];
GetNlist(list, title, bname, !input);
if (input) {
    fprintf(fp, 
     "<caption> All connections between cell types, %d total connections on %d instances."
     "  Cells of type %s have receptors for (at least) <div title=\"%s\">%s</div></caption>\n", 
     total, row_ids.size(), bname, title, list);
   }
else {
    fprintf(fp, 
     "<caption> All connections between cell types, %d total connections on %d instances."
     "  Cells of type %s use neurotransmitter <div title=\"%s\">%s</div></caption>\n", 
     total, row_ids.size(), bname, title, list);
   }
fprintf(fp, 
"<TD> Cell type <TD> %s <TD> Count <TD> Per Inst. <TD> Mean <TD> Conn per Inst<TD> Std. Dev <TD> Dev as %% of Mean<TD> Total %%"
"<TD>7col count<TD>7col mean<TD>7col std<TD>pale count<TD>pale mean<TD>pale std<TD>Yellow count<TD>Yellow mean<TD>Yellow std"
"</TD></TR>\n", input?"Trans":"Receptors");
double sum7 = 0.0;
double weak = 0.0;
for(int ii=0; ii<NumTypes; ii++) {
    int i = sa[ii].index;
    const char *p = CellTypes[i].name;
    if (stats[i].HowMany() > 0) { // if any connections at all, write a row
        double mean, std;
        stats[i].Stats(mean, std);
        fprintf(fp, "<TD BGCOLOR=\"#%06x\"><a href=\"%s-%s.html\"><pre>%s</pre></a>\n", color(p), suffix, p, p);
        char list[128], title[256];
        GetNlist(list, title, p, input);
        fprintf(fp,"<TD><DIV title=\"%s\">%s</DIV>\n", title, list);
        double cpi = mean*double(stats[i].HowMany())/row_ids.size(); // connections per instance
        int color = (cpi > 4.999) ? 0xAAFFAA : (cpi < 0.999 ? 0xFFFFAA : 0xFFFFFF);
        fprintf(fp, "<TD> %d <TD> %.2f <TD> %.2f <TD BGCOLOR=\"#%06x\"> %.2f\n", 
         stats[i].HowMany(), double(stats[i].HowMany())/row_ids.size(), mean, color, cpi );
        if (stats[i].HowMany() == 1)
	    fprintf(fp,"<TD> <TD>");  // no std if only one.
        else {
            fprintf(fp, "<TD> %.2f\n",   std); //std
            fprintf(fp, "<TD> %.2f%%\n", std/mean*100.0);  //frac
            }
        fprintf(fp, "<TD> %.2f%%\n", mean*stats[i].HowMany()/double(total)*100.0);
        PrintEntries(fp, seven [i], &sum7, &weak);
        PrintEntries(fp, pale  [i]);
        PrintEntries(fp, yellow[i]);
        fprintf(fp, "</TD> </TR>\n");
	}
    };
fprintf(fp, "</TABLE>\n");
int ty = TypeIndex(names[row_ids[0]]);
if (CellTypes[ty].sc_index >= 0 && sum7 >= 100)
     fprintf(fp, "Weak tail %.2f of %.2f = %.2f%%\n", weak, sum7, weak/sum7*100.0);
fprintf(fp,"</p>\n");

// Yet another table, optimized for looking at the long tail
vector<vector<string> > cts; // the counts specified as a comma separated string.
vector<vector<int> > totals; // total of above counts
vector<vector<int> > tcount; // How many contribute
vector<int> sum_by_row;
vector<sort_down> any(NumTypes);
int grand_total = 0;
for(int i=0; i<NumTypes; i++)
    any[i].index = i;
for(int i=0; i<row_ids.size(); i++) {
    vector<string> vs(NumTypes);
    vector<int> total_by_type(NumTypes,0);
    vector<int> count_by_type(NumTypes,0);
    int sum = 0;
    for(int j=0; j<rows[i].size(); j++) {
        int sc = dat[i][j];
        if (sc > 0) { // append to string
	    char buf[128];
            sprintf(buf,",%d", sc);
	    int ty = TypeIndex(names[rows[i][j]]);
	    vs[ty] += buf;
	    any[ty].val += dat[i][j];
            grand_total += dat[i][j];
            total_by_type[ty] += sc;
            count_by_type[ty]++;
            sum += sc;
	    }
	}
    cts.push_back(vs);
    totals.push_back(total_by_type);
    tcount.push_back(count_by_type);
    sum_by_row.push_back(sum);
    }
// sort in decreasing order of total connection strength
sort(any.begin(), any.end());
// If "Unknown" is not the last, Move it to last
for(int j=0; j<any.size()-1; j++) {
    if (any[j].index == NumTypes-1) {  // marked unknown
        sort_down save = any[j];
        for(int k=j+1; k<any.size(); k++)
	    any[k-1] = any[k];
        any[any.size()-1] = save;
	break;
	}
    }

fprintf(fp,"<p>\n");
fprintf(fp, "<TABLE BORDER=4 CELLSPACING=4 CELLPADDING=4>\n");
fprintf(fp, "<TR><TD>");
for(int i=0; i<NumTypes; i++) {
    if(any[i].val > 0) {
        int k = any[i].index;
	char linkname[128];
        sprintf(linkname, "%s-%s.html", input ? "out" : "in", CellTypes[k].name);
        int pct = int(any[i].val/grand_total*100.0 + 0.5);
        fprintf(fp, "<TD BGCOLOR=\"#%06x\"><a href=\"%s\">%s(%d%%)</a>\n",
         color(CellTypes[k].name), linkname, CellTypes[k].name, pct);
        }
    }
// write data for Kenyon Cell histograms here
FILE *fkc = fopen("KCdata.txt", "a");
if (fkc == NULL) {
    printf("Could not open KCdata.txt\n");
    exit(42);
    }
fprintf(fp, "</TD></TR>\n");
// Now write each row
for(int i=0; i<row_ids.size(); i++) {
    int color_row = incomplete[i] ? (CompVol.HowMany() > 0 && volumes[i] > 0.9 * CompVol.Mean() ? 0xFFFF99 : 0xFFFFFF) : 0x99FF99;
    fprintf(fp, "<TR><TD BGCOLOR=\"#%06x\">%s", color_row, names[row_ids[i]]);
    fprintf(fkc,"%s %s ", suffix, names[row_ids[i]]);
    for(int j=0; j<NumTypes; j++) {
        if(any[j].val > 0) {
            int k = any[j].index;
            cts[i][k].erase(0,1);  // remove the unwanted initial comma
            int pct = int(double(totals[i][k])/sum_by_row[i]*100.0 + 0.5);  // 0.5 for rounding
            if (totals[i][k] > 0) {
	        fprintf(fp, "<TD>%s (#%d,%d,%d%%)", cts[i][k].c_str(), tcount[i][k], totals[i][k], pct);
                fprintf(fkc,"%s %s ", CellTypes[k].name, cts[i][k].c_str() );
                }
            else
	        fprintf(fp, "<TD> ");
            }
	}
    fprintf(fp, "</TD></TR>\n");
    fprintf(fkc, "\n");
    }
fclose(fkc);

fprintf(fp, "</TABLE>\n");
fprintf(fp,"</p>\n");
fprintf(fp, "</body>\n");
fprintf(fp, "</html>\n");
fclose(fp);

printf("RF: ------------------ %s ------------------\n", names[row_ids[0]] );
// Using the previous table, write a set of commands for computing receptive fields
for(int i=0; i<row_ids.size(); i++) {
    printf("AVG %s", NoSpace(names[row_ids[i]]));
    // Now for the 5 most numerous input types, print contribution
    for(int k=0; k<5; k++) {
	int ty_look_for = any[k].index;
        int sum = 0;
        for(int j=0; j < rows[i].size(); j++) {
	    int ty = TypeIndex(names[rows[i][j]]);
            if (ty == ty_look_for)
		sum += dat[i][j];
	    }
	printf(" %d %s-from-%s", sum, NoSpace(names[row_ids[i]]), CellTypes[ty_look_for].name);
	}
    printf("\n");
    // now figure out how each of them got here
    for(int k=0; k<5; k++) {
	int ty_look_for = any[k].index;
        printf(" AVG %s-from-%s", NoSpace(names[row_ids[i]]), CellTypes[ty_look_for].name);
        for(int j=0; j < rows[i].size(); j++) {
	    int ty = TypeIndex(names[rows[i][j]]);
            if (ty == ty_look_for)
		printf(" %d %s", dat[i][j], NoSpace(names[rows[i][j]]) );
	    }
	printf("\n");
	}
    }

// Now an extra test looking for a cell with a strong connection that is not reported in another cell marked as complete
const int STRONG = 5;
for(int i=0; i<row_ids.size(); i++) {
    for(int j=0; j<rows[i].size(); j++) {
	if (dat[i][j] >= STRONG) { // we found a strong connection.
	    // search for another cell that is (a) complete, and (b) has no corresponding connection
	    int t1 = TypeIndex(names[rows[i][j]]);
	    for(int ii = 0; ii < row_ids.size(); ii++) {
		if (ii == i || incomplete[ii])
		    continue;
		// OK, cell is thought to be complete.  Does it have the same type of connection?
		bool got = false;
		for(int jj = 0; jj < rows[ii].size(); jj++) {
		    if (TypeIndex(names[rows[ii][jj]]) == t1 && dat[ii][jj] > 1)
			got = true;
		    }
		if (!got) {
		    printf("### Very odd.  %d synapses, Cell %s %s %s, but nothing in complete cell %s\n",
		    dat[i][j], names[row_ids[i]], input ? "inputs from" : "outputs to", names[rows[i][j]], names[row_ids[ii]] );
		    }
		}
	    }
	}
    }

}
//---------------------------------- Print out a stereotypy table.
void PrintStereo(FILE *fd, vector<CellStat> &stats, vector<sort_down> &sa, const char *label) {
// Contents go here
fprintf(fd, "<p>%s\n<TABLE BORDER=4 CELLSPACING=4 CELLPADDING=4>\n", label);
fprintf(fd, "<TR><TD> &sigma; inputs<TD>&sigma; outputs<TD>&sigma; volume<TD>&sigma; area<TD>Used/Total<TD>Name</TD></TR>\n");
for(int i=sa.size()-1; i >= 0; i--) { // since array sorts backwards
    int j = sa[i].index;
    // Only print the ones with at least 2 complete cells
    if (CellStats[j].used_stats >= 2) {
        fprintf(fd, "<TR><TD>%.2f%% <TD>%.2f%%<TD> %.2f%%<TD> %.2f%%<TD> %d/%d<TD BGCOLOR=\"#%06x\"> %s</TD></TR>\n", 
         CellStats[j].in_dev*100, CellStats[j].out_dev*100, CellStats[j].volume_dev*100, CellStats[j].area_dev*100, 
         CellStats[j].used_stats, CellStats[j].how_many, CellStats[j].color, CellStats[j].name);
	}
   }
fprintf(fd, "</TABLE>\n</p>\n");
}

// breaks a list up into comma separated parts.  Allocates results so the argument does not need to be saved.
vector<const char *> CommaList(char *list) {
vector<const char *> rslt;
for(char *p = strtok(list, ","); p != NULL; p = strtok(NULL, ",") )
    rslt.push_back(strdup(p));
return rslt;
}

double Inner(vector<double> &a, vector<double> &b) {
if (a.size() != b.size() ){
    printf("Bogus sizes\n");
    exit(42);
    }
double aa = 0.0, bb = 0.0, ab = 0.0;
for(int i=0; i<a.size(); i++) {
    aa += a[i]*a[i];
    bb += b[i]*b[i];
    ab += a[i]*b[i];
    }
//printf("aa = %f bb = %f ab = %f\n", aa, bb, ab);
if (aa*bb <= 0.0)  // should never happen
    return 0.0;
return ab/sqrt(aa*bb);
}

// Figure out odds that the cell with measurements InProfile and OutProfile could be a fragment of a neuron of
// type k.
//
double CouldBeOdds(int k, vector<int> &InProfile, vector<int> &OutProfile, vector<vector<int> > &Strongest)
{
double rslt = 1.0;
// compare the observed outs to the strongest ever seen for that type
for(int j=0; j<NumTypes-1; j++) {
    if (OutProfile[j] <= Strongest[k][j])  //completely consistent with evidence, do nothing
	continue;
    if (Strongest[k][j] == 0) // never seen this before
	rslt *= pow(0.1, OutProfile[j]);  // the stronger the profile, the less likely this could be
    else {                     // we have seen a connection before, but this one is stronger
        double ratio = double(OutProfile[j])/Strongest[k][j];
        rslt *= exp(-ratio);   // the bigger the ratio, the less likely
        }
    }
// do the same for inputs
for(int j=0; j<NumTypes-1; j++) {
    if (InProfile[j] <= Strongest[j][k])  //completely consistent with evidence, do nothing
	continue;
    if (Strongest[j][k] == 0) // never seen this before
	rslt *= pow(0.1, InProfile[j]);  // the stronger the profile, the less likely this could be
    else {                     // we have seen a connection before, but this one is stronger
        double ratio = double(InProfile[j])/Strongest[j][k];
        rslt *= exp(-ratio);   // the bigger the ratio, the less likely
        }
    }
return rslt;
}

// Adds to a strength histogram.  Axis is log, but area added is proportional to strength
void AddToStrengthHisto(vector<double> &histo, double str)
{
// calculate log, base 1.2
double lb12 = log(str)/log(1.2);
// find the int below
double f = floor(lb12);
double alpha = lb12 - f;
int n = int(f);
if (n < 0 || n >= histo.size()-1) {
    printf("very odd in histo %d\n", n);
    exit(42);
    }
// use of alpha made graph smoother, but more confusing
//histo[n]   += (1.0-alpha) * str;
//histo[n+1] += (    alpha) * str;

histo[n] += str;
}

// finds the most common value in an array using a very stupid algorithm.
int FindMode(vector<int> &a, int *HowMany = NULL) {
int mode = 0;
int how_many = 0;
for(int i=0; i<a.size(); i++) {
    // how many equal to entry [i]?
    int N = 0;
    for(int j=0; j<a.size(); j++)
        N += int(a[j] == a[i]);
    if (N > how_many) {
	how_many = N;
        mode = a[i];
	}
    }
if (HowMany != NULL)         // return how many had the mode, if wanted
    *HowMany = how_many;
return mode;
}

void LookForMissing(vector<vector<unsigned char> > &from_to, vector<char *>& names, map<int,int> &rev_map,
 int from, int to, double asnze)
{
// first make sure the cells exist.  One source of a non-existing connection can be a missing cell.
if (from < 0 || to < 0)
    return;
printf("MISSING connection between %s and %s, strength %.2f\n", names[from], names[to], asnze);
// First look for an unnamed fragment that is pre-synaptic to 'to'.  We will see if any of these are adjacent to 'from'.
int body_id_from = rev_map[from];
int body_id_to   = rev_map[  to];
printf("{\"Id1\":%d,\"Names\":\"%s-a-(?)->%s\", \"Id2\":[", body_id_from, names[from], names[to]);
bool first = true;
for(int k=0; k<names.size(); k++) {
    if (from_to[k][to] > 0 && isdigit(names[k][0]) ){
	printf("%c%d", first? ' ' : ',', rev_map[k]);
        first = false;
	}
    }
printf("]}\n");

// now look for fragments post-synaptic to 'from'.
printf("{\"Id1\":%d,\"Names\":\"%s->()-a-%s\", \"Id2\":[", body_id_to, names[from], names[to]);
first = true;
for(int k=0; k<names.size(); k++) {
    if (from_to[from][k] > 0 && isdigit(names[k][0]) ){
	printf("%c%d", first? ' ' : ',', rev_map[k]);
        first = false;
	}
    }
printf("]}\n");
}

double sqr(double s){ return s*s;}

// Process and write a file for synapse overlaps.  If M > 0, it's for a specific medulla layer.  If M < 0, it's
// for all layers.
void ProcessSynOverlaps(const char *from, const char *to, double *synapse, double *overlap, FILE *overview, int m)
{
char filename[256];
if (m > 0)
    sprintf(filename, "scatter%d/%s-%s", m,  from, to);
else
    sprintf(filename, "scatter/%s-%s", from, to);
FILE *fp = fopen(filename, "w");
if (fp == NULL) {
    printf("Could not %s for write\n", filename);
    exit(42);
    }
MeanStd s,o;
MeanStd d;  // density
for(int i=0; i<7; i++) {
    fprintf(fp, "%.2f %.2f\n", overlap[i], synapse[i]);
    s.Element(synapse[i]);
    o.Element(overlap[i]);
    if (overlap[i] > 0)
        d.Element(synapse[i]/overlap[i]);
    }
fprintf(fp, "\n");
fclose(fp);

// calculate the coefficent of correlation, r, between area and synapse count.
double sum1 = 0.0;
double vs = 0.0, va = 0.0;
for(int i=0; i<7; i++) {
    sum1 += (synapse[i]-s.Mean()) * (overlap[i] - o.Mean());
    vs   += sqr(synapse[i] - s.Mean());
    va   += sqr(overlap[i] - o.Mean());
    }
double r;
if (va > 0 && vs > 0)
    r = sum1/sqrt(vs)/sqrt(va);
else {
    printf("One coord had no variance? %f %f %s %s\n", va, vs, from, to);
    r = 0.0;
    }

double density = (d.HowMany() > 0) ? d.Mean() : 1.0;
fprintf(overview, "%.2f %.2f %.2f %.3f \"%s-%s\" %.3f \n", o.Mean(), s.Mean(), s.Std(), density, from, to, r);

// Find best fit to synapse count as a function of area
vector<double> area(7), syn(7);
double min_area =  1e30;
double max_area = -1e30;
for(int i=0; i<7; i++) {
    area[i] = overlap[i];
    syn [i] = synapse[i];
    min_area = min(min_area, area[i]);
    max_area = max(max_area, area[i]);
    }
Fitab fit(area, syn);
printf("  %s fit %s, layer %d: Syn as a function of area %f %f\n", from, to, m, fit.a, fit.b);
double y = fit.a +fit.b * min_area;
printf("%.3f %.3f \"\"\n", min_area, y);
y =        fit.a +fit.b * max_area;
printf("%.3f %.3f \"%s-%s\"\n\n", max_area, y, from, to);
}

// See if there is some combination of two vectors that has less variance than either alone.
double LookForBetter(int i, int j, vector<double> &v1, vector<double> &v2, vector<const char *> &names, 
vector<sort_down> &rslts, char op, int pass, bool print_recip = false)
{
double amt = 100.0;  // amount (if any) by which it is better than either of its two parents
                    // Return a value of 100 for did not compute
int N = names.size();
int nz1 = 0;  // non-zero counts
int nz2 = 0;
MeanStd mv1, mv2;
for(int k=0; k<7; k++) {
    mv1.Element(v1[k]);
    mv2.Element(v2[k]);
    if (v1[k] != 0 ) nz1++;
    if (v2[k] != 0 ) nz2++;
    }
if (nz1 == 0 || nz2 == 0)               // both directions must have at least one non-zero
    return amt;                         // or else return the 'bogus' flag
double cn = 0.0, cd1 = 0.0, cd2 = 0.0;  // for correlation
for(int k=0; k<7; k++) {
    double d1 = v1[k] - mv1.Mean();
    double d2 = v2[k] - mv2.Mean();
    cn += d1*d2;
    cd1 += d1*d1;
    cd2 += d2*d2;
    }
double corr = (abs(cd1*cd2) <1e-6) ? 0.0 : cn/sqrt(cd1*cd2);  // in case get all the same integer in
if (pass == 2)
    printf("Correlation is %.4f\n", corr);
double fm1 = mv1.Std()/mv1.Mean();
double fm2 = mv2.Std()/mv2.Mean();
if (nz1 == 7 && nz2 == 7) {   // try the ratio
    MeanStd m;
    vector<double> ratios(7);
    for(int k=0; k<7; k++) {
        double r;
        if (op == '+')
	    r = v1[k] + v2[k];
        else if (op == '/')
	    r = v1[k] / v2[k];
        else {
	    printf("Bad operator '%c'\n", op);
	    exit(42);
	    }
	m.Element(r);
	ratios[k] = r;
	}
    double fm = m.Std()/m.Mean();
    // It is possible that fm1 or fm2 is 0, if we were given a vector of all the same integer.
    if (fm1 > 1e-6 && fm2 >= 1e-6)
        amt = fm/min(fm1, fm2);
    else
        amt = 100.0;  // for a bad result.
    if (print_recip && (i/N) == (j%N) && (i%N) == (j/N)) { // it's reciprocal.
	printf("\nRecip: Operator %c compared to either for %s->%s / %s->%s by %.3f .  Corr=%.4f:\n", 
         op, names[i/N], names[i%N], names[j/N], names[j%N], amt, corr);
	for(int k=0; k<7; k++)
	    printf("     %6.1f %c %6.1f  = %6.3f\n", v1[k], op, v2[k], ratios[k]);
        printf("     ------   ------    ------\n");
	printf("mean %6.1f   %6.1f    %6.3f\n", mv1.Mean(), mv2.Mean(), m.Mean() );
	printf(" std %6.2f   %6.2f    %6.3f\n", mv1.Std(),  mv2.Std(),  m.Std()  );
	printf("norm %6.3f   %6.3f    %6.3f\n", fm1,        fm2,        fm       );
        }
    // For trying direct correlation rather than our dev/min of input dev metric.
    // Will return 0 (perfectly correlated) to 20 (perfectly anti-correlated).
    // We want 1.0 to be a good result, going down to 0 as a perfect result
    bool use_corr = true;  // instead of ratio of normalized deviation
    if (use_corr)
        amt = op == '/' ? (1.0-corr)*10.0 : (1.0+corr)*10.0;  // for '+', want anti-correlated
    if (amt < 1) {
        if (pass == 1)   // in pass 1, record it
            rslts.push_back(sort_down(i,j,amt));
        else {           // in pass 2, write it out
	    printf("Operator %c better than either for %s->%s / %s->%s by %.3f:\n", op, names[i/N], names[i%N], names[j/N], names[j%N], amt);
	    for(int k=0; k<7; k++)
	        printf("%6.1f %6.1f %6.3f\n", v1[k], v2[k], ratios[k]);
            }
	}
    }
return amt;
}

// Random numbers, Numerical recipes with superficial changes
#define SQR(u) ((u)*(u))

struct Ranq1 {
    unsigned long long v;
    Ranq1(unsigned long long j): v(4101842887655102017LL) {
	v ^= j;
	v = int64();
	}
    inline unsigned long long int64() {
	v ^= v >> 21; v ^= v << 35; v ^= v >> 4;
	return v * 2685821657736338717LL;
	}
    inline double doub(){return 5.42101086242752217E-20 * int64();}
    };

struct NormalDev : Ranq1 {
    double mu, sig;
    NormalDev(double mmu, double ssig, unsigned long long i) :
	Ranq1(i), mu(mmu), sig(ssig){}
    double dev() {
	double u,v,x,y,q;
        do {
	    u = doub();
	    v = 1.7156*(doub() - 0.5);
            x = u - 0.449871;
            y = fabs(v) + 0.386595;
            q = SQR(x) + y*(0.19600*y-0.25472*x);
	} while (q > 0.27597 && (q > 0.27846 || SQR(v) > -4.*log(u)*SQR(u)));
	return mu + sig*v/u;
	}
    };

class tail {
  public:
    short dat[10];
    tail() { for(int i=0; i<10; i++) dat[i] = 0;}
    void add(double d) { if (d < 1.0){ int i = int(d*10.0); dat[i]++;};}
    };

// Look for correlations between the columns
void SearchForCorrelation(vector<vector<double> > dd, int Nrac, vector<const char *> &all_col_names, char op, const char *fname)
{

vector<sort_down> better;

// With so many possibilities, we are sure to get some false positives.
//So make a fake version for testing the null hypothesis.
//Actually, make 100, or 10000, versions, and keep track of their mean and standard deviation.
//We can use this to see how unlikely any real results could be
vector<int> envelope(101,0);
vector< vector< MeanStd> > FakeStats(Nrac*Nrac, vector< MeanStd>(Nrac*Nrac) );
vector< vector< tail> >    FakeTails(Nrac*Nrac, vector<    tail>(Nrac*Nrac) );
MeanStd StatsOfFake;
MeanStd FakeSameCell;
MeanStd FakeReciprocal;
for(int trial=0; trial<100; trial++) {   // Use 10000 for real resultss
    vector<vector<double> > save_dd = dd;
    for(int i=0; i<Nrac*Nrac; i++) {
	// replace each row with random numbers with same mean and std deviation.
	MeanStd m;
	for(int k=0; k<7; k++)
	    m.Element(dd[i][k]);
	NormalDev norm(m.Mean(), m.Std(), trial * 1000000 + i /* random seed */);
	if (m.Sum() > 1) {
	    double rnd = 0.5;
	    for(int k=0; k<7; k++)
		dd[i][k] = int(max(norm.dev(), 0.0) + rnd); // negative numbers cause chaos.  Delete them even if theoretically dubious.
	    }
	}
    better.clear();
    vector<int>ratio_hist2(101,0);
    int Nsamecell = 0;
    int Nrecip   = 0;
    for(int j=0; j<Nrac*Nrac; j++)
	for(int k=0; k<Nrac*Nrac; k++)
	    if (j != k) {
		double amt = LookForBetter(j, k, dd[j], dd[k], all_col_names, better, op, 1 /* pass 1, just record */);
		int f1 = (j/Nrac); // from and to of both
		int t1 = (j%Nrac);
		int f2 = (k/Nrac);
		int t2 = (k%Nrac);
		bool recip = (f1 == t2 && f2 == t1);
		bool same_cell = (f1 == t2 || f1 == f2 || t1 == f2 || t1 == t2);
		if (amt < 1.0 && recip)
		    Nrecip++;
		if (amt < 1.0 && same_cell)
		    Nsamecell++;
                if (amt != 100.0) {
                    FakeStats[j][k].Element(amt);
                    FakeTails[j][k].add(amt);     // add to histogram
		    }
		}
    printf("Null hypothesis (fake data, all pairs) gives %d results with score < 1.\n", better.size() );
    StatsOfFake.Element(better.size());
    FakeSameCell.Element(Nsamecell);
    FakeReciprocal.Element(Nrecip);
    for(int i=better.size()-1; i>=0; i--)
	ratio_hist2[int(100.0*better[i].val)]++;
    for(int i=0; i<=100; i++)
	envelope[i] = max(envelope[i], ratio_hist2[i]);
    dd = save_dd;
    }
printf("From %d fake example, on average %.3f (std %.3f) has scores of less than 1.\n", StatsOfFake.HowMany(), StatsOfFake.Mean(), StatsOfFake.Std() );
printf("From %d fake example, on average %.3f (std %.3f) involved a common cell.\n", FakeSameCell.HowMany(), FakeSameCell.Mean(), FakeSameCell.Std() );
printf("From %d fake example, on average %.3f (std %.3f) were reciprocal.\n", FakeReciprocal.HowMany(), FakeReciprocal.Mean(), FakeReciprocal.Std() );
// Now that we have some idea what to expect, try real data.
better.clear();
MeanStd all_pairs, recip_pairs, same_cell_pairs;
for(int j=0; j<Nrac*Nrac; j++)
    for(int k=0; k<Nrac*Nrac; k++)
	if ((op == '+') ? k<j : k != j ) {  // op '+' needs only k,j since j,k will be identical; op '/' needs all k!=j
	    double amt = LookForBetter(j, k, dd[j], dd[k], all_col_names, better, op, 1 /* pass 1 */, true /* print recip */);
            int f1 = (j/Nrac); // from and to of both
            int t1 = (j%Nrac);
            int f2 = (k/Nrac);
            int t2 = (k%Nrac);
            bool recip = (f1 == t2 && f2 == t1);
            bool same_cell = (f1 == t2 || f1 == f2 || t1 == f2 || t1 == t2);
            if (amt < 100 && recip)
		recip_pairs.Element(amt);
            if (amt < 100 && same_cell)
		same_cell_pairs.Element(amt);
            if (amt < 100)
		all_pairs.Element(amt);
	    }
printf("Average over all %d pairs is %.3f, over %d same cell pairs %.3f, over %d reciprocal pairs is %.3f\n", 
 all_pairs.HowMany(), all_pairs.Mean(), same_cell_pairs.HowMany(), same_cell_pairs.Mean(), recip_pairs.HowMany(), recip_pairs.Mean() );
// Replay list, this time sorted by metric
sort(better.begin(), better.end());
printf("------------------- Replay list, %d now sorted ------------------------\n", better.size());
vector<int>ratio_hist1(101,0);
for(int i=better.size()-1; i>=0; i--) {
    printf("Expecting %.3f\n", better[i].val);
    int j = better[i].index;
    int k = better[i].index2;
    ratio_hist1[int(100.0*better[i].val)]++;
    double amt = LookForBetter(j, k , dd[j], dd[k], all_col_names, better,  op, 2 /* pass 2 */);
    double me = FakeStats[j][k].Mean();
    double st = FakeStats[j][k].Std();
    double z;  // number of stds out
    if (abs(st) < 1e-8) {
	printf("Std too small %f (mean %.3f N=%d\n", st, me, FakeStats[j][k].HowMany());
        z = 1.0;
	}
    else
        z = (me - amt)/st;
    int N = FakeStats[j][k].HowMany();
    const char *flag =  (N >= 90 && abs(z) >= 3.0) ? "###" : "";
    printf("  Bing! entry %.3f is %.2f standard deviations (each %.3f) from mean %.3f of %d samples %s\n", amt, z, st, me, N, flag);
    double cumul = 0.0;
    for(int m=0; m<10; m++) {
	printf("%5d", FakeTails[j][k].dat[m]);
        double bin_bot =  m   *0.1;
        double bin_top = (m+1)*0.1;
        if (bin_top < amt)
	   cumul += FakeTails[j][k].dat[m];
        if (bin_bot <= amt && amt <= bin_top)
	   cumul += (amt - bin_bot) * 10.0 * FakeTails[j][k].dat[m];
	}
    printf(" :  Cumulative %.2f\n", cumul);
    }

FILE *fp = fopen(fname, "w");
if (fp == NULL) {
    printf("Could not open '%s' for write\n", fname);
    exit(42);
    }
for(int i=0; i<=100; i++) {
     printf(    "%4d %6d %4d\n", i, ratio_hist1[i], envelope[i]);
    fprintf(fp, "%4d %6d %4d\n", i, ratio_hist1[i], envelope[i]);
    }
fclose(fp);
}

// Routine for fitting a plane.  Styled after NR example.
VecDoub Plane2D(VecDoub_I &xx) {
VecDoub ans(3);
double x = xx[0];
double y = xx[1];
ans[0] = 1;
ans[1] = x;
ans[2] = y;
return ans;
}

// Tells if point c is on the left side of the vector a->b by using the cross product
bool LeftSide(const Point &a, const Point &b, const Point &c)
{
return (b.x-a.x)*(c.y-a.y) - (b.y-a.y)*(c.x-a.x) > 0;
}

bool LinesCross(const Point &p1, const Point &p2, const Point &p3, const Point &p4, Point *rslt = NULL)
{
double a = p2.x - p1.x;
double b = p3.x - p4.x;
double c = p2.y - p1.y;
double d = p3.y - p4.y;
double det = a*d - b*c;

double e = p3.x - p1.x;
double f = p3.y - p1.y;
double alpha = double(d*e - b*f) / double(det);
double beta = double((-c)*e + a*f) / double(det);
bool cross = 0.0 < alpha && alpha < 1.0 && 0.0 < beta && beta < 1.0;
if (cross) {
    printf("Aha!  got crossing (%f %f) to (%f %f) crosses (%f %f) to (%f %f), alpha=%f beta=%f\n",
     p1.x, p1.y, p2.x, p2.y, p3.x, p3.y, p4.x, p4.y, alpha, beta);
    double ax = alpha*p2.x + (1-alpha)*p1.x;
    double ay = alpha*p2.y + (1-alpha)*p1.y;
    double bx =  beta*p4.x + (1-beta) *p3.x;
    double by =  beta*p4.y + (1-beta) *p3.y;
    printf("Intersect (%f %f) and (%f %f)\n", ax, ay, bx, by);
    if (rslt != NULL) {
        rslt->x = ax;
        rslt->y = ay;
	}
    }
return cross;
}

// Go through all the tbars, looking for precursor cells.
// // ty = type of T4 cell.  May be "fb" (front-back) or "bf" or "ud (up-down) or "du".
void MapXYs(vector<TBar> &tbs, int what, const char *ty, map<int,int> &name_map, vector<char *> &names)
{
vector<Partner> p;
Point center(0.0,0.0);
int count = 0;
for(int i=0; i<tbs.size(); i++) {
    for(int j=0; j<tbs[i].partners.size(); j++) {
	if (tbs[i].partners[j].body_id == what) {
            Point pt = tbs[i].partners[j].pt;
	    printf("Mike: %6d from %6d at x,y=%.1f %.1f %s\n", what, tbs[i].body_id, pt.x, pt.y,
             names[name_map[tbs[i].body_id]]);
            Partner pp = tbs[i].partners[j];
            pp.body_id = tbs[i].body_id;
            p.push_back(pp);
	    count++;
            center.x += pt.x; center.y += pt.y;
	    }
	}
    }
center.x /= count; center.y /= count;  // this is the overall center

printf("Sorted by source id\n");
sort(p.begin(), p.end() );
// If option is enabled, find the center of the strongest connection instead of overall geometric center
int strongest = 0;  // strength of strongest connection
bool use_strongest = true;
int last = -1;
for(int i=last+1; use_strongest && i<p.size(); i=last+1) {
    for(last=i; last < p.size()-1 && p[last+1].body_id == p[i].body_id; last++)
	;
    int N = last - i + 1;
    if (N > strongest) {
        strongest = N;
        center = Point(0.0, 0.0);
        for(int j=i; j<=last; j++) {
	    center.x += p[j].pt.x; center.y += p[j].pt.y;
	    }
        center.x /= N; center.y /= N;
	}
    }
// OK, move to 0,0
for(int i=0; i<p.size(); i++) {
    p[i].pt.x -= center.x;
    p[i].pt.y -= center.y;
    }

vector <double> SumXbyType(NumTypes, 0.0);
vector <double> SumYbyType(NumTypes, 0.0);
vector <int>    SumNbyType(NumTypes, 0.0);
last = -1;
char fname[128];
sprintf(fname, "T4/%s/other", ty);    // other cell types
FILE *fother = fopen(fname, "w");
if (fother == NULL) {
    printf("Could not open '%s'\n", fname);
    exit(42);
    }

for(int i=last+1; i<p.size(); i=last+1) {
    for(last=i; last < p.size()-1 && p[last+1].body_id == p[i].body_id; last++)
	;
    const char *Type = BaseName(names[name_map[p[i].body_id]]);
    int ti = TypeIndex(Type);
    double sumx = 0.0, sumy = 0.0;
    for(int j=i; j<=last; j++) {
        printf("From %6d at %.1f %.1f\n", p[j].body_id, p[j].pt.x, p[j].pt.y);
	sumx += p[j].pt.x;
        sumy += p[j].pt.y;
	}
    printf("\n");
    SumXbyType[ti] += sumx;
    SumYbyType[ti] += sumy;
    SumNbyType[ti] += (last - i +1);;
    sumx /= (last-i+1);
    sumy /= (last-i+1);
    double r = sqrt(double(last-i+1));  // so area of circle is proportional to weight
    // print the 'spokes' from the center. 
    sprintf(fname, "T4/%s/%s.lines", ty, Type);    // main cell types
    FILE *flines = fopen(fname, "a");
    if (flines == NULL) {
	printf("Could not open '%s'\n", fname);
	exit(42);
	}
    for(int j=i; j <= last; j++) 
        fprintf(flines, "%.1f %.1f\n%.1f %.1f\n\n", p[j].pt.x, p[j].pt.y, sumx, sumy);
    fclose(flines);
    // if a T4 cell, print the first char of direction ( [4] )  at the center.. If it's 't'
    // then the direction is unknown.  Print '?'.
    if (strcmp(Type, "T4") == 0) {
	sprintf(fname, "T4/%s/T4.labels", ty);
	FILE *flabels = fopen(fname, "a");
	if (flabels == NULL) {
	    printf("Could not open '%s'\n", fname);
	    exit(42);
	    }
        char dir = names[name_map[p[i].body_id]][4];
	fprintf(flabels, "%.1f %.1f \"%c\"\n", sumx, sumy, dir == 't' ? '?' : dir);
        fclose(flabels);
	}
    sprintf(fname, "T4/%s/%s", ty, Type);
    FILE *fp = fopen(fname, "a");
    if (fp == NULL) {
	printf("cannot open '%s'\n", fname);
	exit(42);
	}
    for(int j=i; j<= last; j++)
	fprintf(fp,"%.1f %.1f 10.0\n", p[j].pt.x, p[j].pt.y);
    fprintf(fp, "%.1f %.1f %.1f\n", sumx, sumy, 40.0*r);
    fclose(fp);
    }
// write a file for the center of gravity of each type
for(int i=0; i<NumTypes-1; i++) {
    int N = SumNbyType[i];
    if (N <= 0)
	continue;
    sprintf(fname, "T4/%s/%s.labels", ty, CellTypes[i].name);
    FILE *flabels = fopen(fname, "a");
    if (flabels == NULL) {
	printf("Could not open '%s'\n", fname);
	exit(42);
	}
    double x = SumXbyType[i]/N;
    double y = SumYbyType[i]/N;
    fprintf(flabels, "%.1f %.1f \"%s\"\n", x, y, CellTypes[i].name);
    fclose(flabels);
    }
}

void PerpBisect(double x1, double y1, double x2, double y2)
{
double slope = (y2-y1)/(x2-x1);
double bis   = -1.0/slope;
double midx = (x1+x2)/2.0;
double midy = (y1+y2)/2.0;
double len = abs(bis) > 2 ? 150.0 : 500;
// find the x span of a line of length 'len'.
double xspan = sqrt(len*len/(1.0+bis*bis));
double r1x = midx - xspan/2;
double r1y = midy - (xspan/2)*bis;
double r2x = midx + xspan/2;
double r2y = midy + (xspan/2)*bis;
printf("%.1f %.1f\n%.1f %.1f\n\n", r1x, r1y, r2x, r2y);

}

// Trace back inputs starting with Body ID 'root'.
// Connectivity matrix is 'from_to'.
// 'wanted' is the cell types you want to trace.
// 'trace' is false for those cells you do not want to trace back from (typically inputs)
void CircuitBacktrace(
 int root, vector<vector<unsigned char> > &from_to, vector<char *> &names, vector<bool> wanted, vector<bool> trace, int threshold)
{
printf("Michael:\n");
vector<bool> cell_used(names.size(), false);  // to avoid loops
vector<int> look_at;
look_at.push_back(root);  // for now, start with one T4
vector<int> more;

for(int hops=0; look_at.size() > 0; look_at = more) {
    printf("----------- Cells %d hops from target----------\n", hops++);
    more.clear();
    vector<bool> req(names.size(), false);  // requested already on this pass?
    for(int k=0; k<look_at.size(); k++) 
	req[look_at[k]] = true;             // already being considered
    for(int k=0; k<look_at.size(); k++) {
        int n = look_at[k];
        printf(" Cell: %s\n", names[n] );
        cell_used[n] = true;
        // add all inputs from cells of wanted types.
        for(int i=0; i<names.size(); i++) {
            int str = from_to[i][n];
            if (str <= threshold)
		continue;
            int ti = TypeIndex(names[i]);
	    if(wanted[ti]) {
		printf("  Input from %-12s, %3d synapses\n", names[i], str );
		// if we want to trace, and not requested already, add to list    
		if (trace[ti] && (!cell_used[i]) && !req[i]) {
		    more.push_back(i);
		    req[i] = true;  // now been requested
		    }
		}
	    }
	}
    }
printf("------------------------ End of Circuit --------------------------\n");
}
VecDoub FitALinear(const Doub x) {
VecDoub ans(2);
ans[0] = 1.0;
for(int i=1; i<2; i++)
    ans[i] = ans[i-1] * x;
return ans;
}

// Compute distance between partners
double PartnerDist(Partner &a, Partner &b)
{
double dx = a.pt.x - b.pt.x;
double dy = a.pt.y - b.pt.y;
double dz = a.z - b.z;
return sqrt(dx*dx + dy*dy + dz*dz);
}

bool KeepSynapse(bool convergent)
{
return true;
}

vector<Point3d> ReadPajekFile(const char *name) {
vector<Point3d> rslt;
FILE *fp = fopen(name,"r");
if (fp == NULL) {
    printf("Could not open '%s'\n", name);
    exit(42);
    }
char junk[128];
int N;
fscanf(fp, "%s %d", junk, &N);
printf("Opened '%s', got %d entries\n", name, N);
for(int i=0; i<N; i++) {
    double x,y,z;
    int j;
    fscanf(fp, "%d %s %lf %lf %lf", &j, junk, &x, &y, &z);
    rslt.push_back(Point3d(x,y,z));
    }
fclose(fp);
return rslt;
}

// Finds the closest point (that is labelled), then returns 4th char (so KC-s returns 's'),
// for example.  If the next closest is the same, upper case it.
char FindMembership(int k,  vector<Point3d> &rslts, vector<int> &KCs, vector<char *> &names) {
char rslt = 'X';
vector<sort_down> s;
for(int i=0; i<rslts.size(); i++) {
    if (strcmp(names[KCs[i]], "KC-any") == 0)
	continue;
    double dx = rslts[i].x - rslts[k].x;
    double dy = rslts[i].y - rslts[k].y;
    double dz = rslts[i].z - rslts[k].z;
    double extra = (i == k) ? 0.0 : 1.0e-6;  // make sure self-reference sorts first
    s.push_back(sort_down(i, sqrt(dx*dx + dy*dy + dz*dz)+extra));
    };
sort   (s.begin(), s.end());
reverse(s.begin(), s.end());
if (strcmp(names[KCs[k]], "KC-any") != 0 && s[0].index != k)
    printf("How can this happen?? Closest is not itself? k=%d [0].index=%d dist=%f\n", k, s[0].index, s[0].val);
//printf("Closest to %d are:", k);
//for(int j=0; j<4; j++)
    //printf(" %s", names[KCs[s[j].index]] );
//printf("\n");
rslt = names[KCs[s[0].index]][3];
// If the second closest has the same name, upper case (more EMPHATIC)
if (strcmp(names[KCs[s[0].index]], names[KCs[s[1].index]]) == 0)
    rslt = toupper(rslt);
return rslt;
}

void FindCommonInputs(const char *name, vector<int> &KCs, vector<char *>&names, 
 vector<vector<unsigned char> > &from_to, map<int, int> &rev_map) {
vector<int> outs;
for(int i=0; i<names.size(); i++)
    if (strncmp(names[i],name,strlen(name)) == 0) {
	printf("Neuron %s counts for %d and %.4f\n", names[i], (1 << outs.size()), pow(0.01, outs.size()+1) );
	outs.push_back(i);
	}
printf("found %d of %s\n", outs.size(), name );
// now, for every KC, find all the outs it connects to, and generate a bitwise vector
vector<sort_down> sa;
for(int i=0; i<KCs.size(); i++) {
    double sum = 0;
    for(int j=0; j<outs.size(); j++) {
	int str = from_to[KCs[i]][outs[j]];
	if (str > 0)
	    sum += (1 << j) + str * pow(0.01, (j+1));
	}
    sa.push_back(sort_down(i, sum));
    }
sort(sa.begin(), sa.end() );
for(int i=0; i<sa.size(); i++) {
    if (i > 0 && int(sa[i].val) != int(sa[i-1].val))
	printf("\n");
    printf("Shinya: %4d %9d %12s %.4f\n", i, rev_map[KCs[sa[i].index]], names[KCs[sa[i].index]], sa[i].val );
    }
}

int main(int argc, char **argv)
{

printf("Starting!\n");
Json::Value syn_root;   // will contain the root value after parsing synapse file
Json::Value bdy_root;   // will contain the root value after parsing    body file
Json::Value book_root;  // will contain the root value after parsing the bookmarks file
Json::Value small_body_root;   // will contains the root value after parsing  small_body file
Json::Value side_body_root;   // will contains the root value after parsing side_body file
Json::Reader reader;
const int N=90000000;    // create buffer
vector<char> huge(N, 0);

//-------------------------  Read and parse the body JSON file -------------------------------------------
char *fname = strdup("annotations-body.json");
printf("Starting to read file %s\n", fname); fflush(stdout);
FILE *fd = fopen(fname, "r");
if (fd == NULL) {
    printf("Cannot open '%s' for read.\n", fname);
    return 42;
    }
if(fread(&huge[0], 1, N, fd) == N) {
    printf("File too big!\n");
    return 42;
    }
fclose(fd);
printf("JSON body file read, starting to parse.\n"); fflush(stdout);
bool parsingSuccessful = reader.parse(&huge[0], bdy_root );
if ( !parsingSuccessful )
{
    // report to the user the failure and their locations in the document.
    std::cout  << "Failed to parse configuration\n"
               << reader.getFormatedErrorMessages();
    exit(42);
    }

//----------------------------- read the synapse file ------------------------------------------
//
// in this first pass, we just want to make sure we have all the body IDs we will need
map<int,int> name_map;  // map of body ID to name
map<int,int>  rev_map;  // and the reverse
vector<char*> names;    // These are sequential.

fname = strdup("synapse.json");
printf("Starting to read file %s\n", fname); fflush(stdout);
fd = fopen(fname, "r");
if (fd == NULL) {
    printf("Cannot open '%s' for read.\n", fname);
    return 42;
    }
for(int i=0; i<huge.size(); i++)
    huge[i] = '\0';
if(fread(&huge[0], 1, N, fd) == N) {
    printf("File too big!\n");
    return 42;
    }
fclose(fd);
printf("Synapse annotation file read, starting to parse.\n"); fflush(stdout);
parsingSuccessful = reader.parse(&huge[0], syn_root );
if ( !parsingSuccessful )
{
    // report to the user the failure and their locations in the document.
    std::cout  << "Failed to parse configuration\n"
               << reader.getFormatedErrorMessages();
    exit(42);
    }

printf("parsed file...\n");
Json::Value md = syn_root["metadata"];
std::string des = md["description"].asString();
int         ver = md["version"].asInt();

printf("Description: '%s', version %d\n", des.c_str(), ver );

// OK, read and parsed both files.  Start by going through BODY file.  Need to consolidate all CT1
// neurons.  This is a special case where they look like separate body IDs, but we know (from external
// evidence) they are all the same neuron, one that joins outside our valume.  Create a set of ints,
// CT1s, of all the CT1 body IDS.

set<int> CT1s;
bool combine_CT1s = false;
Json::Value data = bdy_root["data"];
if (combine_CT1s) {
for(unsigned int i=0; i<data.size(); i++) {
    //printf("getting entry %d\n", i);
    Json::Value v = data[i];
    // read the different things that can be in 'data'.  For now only "body ID" and "name"
    int t = v["body ID"].asInt();
    printf("Body ID %d\n", t);
    Json::Value part = v["name"];
    if (!part.isNull()) {
	string name = part.asString();
	printf("Name %s\n", name.c_str() );
        if (strncmp(name.c_str(), "CT1", 3) == 0) {
	    printf("Found %s, ID %d\n", name.c_str(), t);
            CT1s.insert(t);
	    }
	}
    }
}

// Go through the synapse file, extracting all body IDs.
data = syn_root["data"];
int Tcount = 0;
int Scount = 0;
for(unsigned int i=0; i<data.size(); i++) {
    //printf("getting entry %d\n", i);
    Json::Value v = data[i];
    // read the different things that can be in 'data'.  For now only {"T-bar","partner"} is known
    Json::Value t = v["T-bar"];
    Json::Value part = v["partners"];
    if (!t.isNull()) {
        Json::Value lo = t["location"];
        TBar tb;
	tb.body_id     = t["body ID"].asInt();
        if (CT1s.find(tb.body_id) != CT1s.end())   // If any of the CT1s, replace by first
	    tb.body_id = *(CT1s.begin());
        if (name_map.find(tb.body_id) == name_map.end()) {
	    name_map[tb.body_id] = names.size();
            rev_map[names.size()] = tb.body_id;
            names.push_back(NULL);
	    }
        Tcount++;
	for(int j=0; j<part.size(); j++) {
	    Partner prt;
	    prt.body_id    = part[j]["body ID"].asInt();
	    if (CT1s.find(prt.body_id) != CT1s.end())   // If any of the CT1s, replace by first
		prt.body_id = *(CT1s.begin());
            if (name_map.find(prt.body_id) == name_map.end()) {
		name_map[prt.body_id] = names.size();
                rev_map[names.size()] = prt.body_id;
		names.push_back(NULL);
	        }
            Scount++;
	    }
	}
    }
printf("Found %d unique IDs, %d T-bars, %d PSDs\n", name_map.size(), Tcount, Scount);

//---------------------- now go through the JSON body file again
vector<char> leaves(name_map.size(), 'U');  // 'U' for unknown, 't' for tangential, 'o' for orphan
fd = fopen("AllNamedBodies.json","w");  // write a file of scripts to get volumes and areas
if (fd == NULL) {
    printf("Could not open script file 'AllNamedBodies.json'\n");
    return 42;
    }
fprintf(fd, "{\"data\":[\n");
int NPrint = 0;
int Mi1 = -1;
int Ncount = 0;
data = bdy_root["data"];
for(unsigned int i=0; i<data.size(); i++) {
    //printf("getting entry %d\n", i);
    Json::Value v = data[i];
    // read the different things that can be in 'data'.  For now only "body ID" and "name"
    int t = v["body ID"].asInt();
    printf("Body ID %d\n", t);
    // Check to make sure it does not exist, then add it to the map
    bool in_the_map = name_map.find(t) != name_map.end();
    //if (!in_the_map)
	//printf("ID %d has no synapses??\n", t);
    Json::Value part = v["name"];
    if (!part.isNull()) {
	string name = part.asString();
        printf("Body name %s\n", name.c_str() );
        if (!in_the_map) {
	    printf("Named body %s has no synapses??\n", name.c_str() );
	    continue;
	    }
        if (name == "Mi1 home")
	    Mi1 = name_map[t];
        names[name_map[t]] = strdup(name.c_str() );
        //printf("echo looking for '%s'\n", name.c_str() );
        if (NPrint++ != 0)
	    fprintf(fd, ",");
        fprintf(fd, "%d", t);
	if ((NPrint % 10) == 0)
	    fprintf(fd, "\n");
        Ncount++;
	}
    part = v["comment"];
    if (!part.isNull()) {
        string val = part.asString();
	if (strstr(val.c_str(), "tangent") != NULL)
	    leaves[name_map[t]] = 't';
	if (strstr(val.c_str(), "orphan" ) != NULL)
	    leaves[name_map[t]] = 'o';
	}
    }
fprintf(fd, "]}\n");
fclose(fd);
printf("Found %d unique IDs, %d T-bars, %d PSDs, %d named bodies\n", name_map.size(), Tcount, Scount, Ncount);
// -------------------- Read the bookmark file and generate column and plane info-----------
vector<vector<vector<double> > > Centers(11, vector<vector<double> >(7, vector<double>(2)));
vector<vector<double> > ZofCenters(11, vector<double>(7));
// ------------------------ Read info about small bodies -------------------------------
// ----------------------- Assign names to all un-named --------------------------------------
map<int,int>::iterator it;
for(it = name_map.begin(); it != name_map.end(); it++) {
    if (names[it->second] == NULL) {
	char buf[128];
        if (leaves[it->second] != 'U')  // if status is known, append it
	    sprintf(buf,"%d%c", it->first, leaves[it->second] );
        else                            // otherwise, just set it equal to the number
            sprintf(buf,"%d", it->first);
        names[it->second] = strdup(buf);
	}
    }
// ------------------------Read the file of sizes and surface area ------------------------------
vector<double> areas(names.size(), -1.0);
vector<double> volumes(names.size(), -1.0);
fd = fopen("AllStats.txt", "r");
if (fd == NULL) {
    printf("No size file present\n");
    }
else {
    // get all the triples.
    int N = -1;
    fscanf(fd, "%d", &N);
    printf("Reading %d sizes and areas\n", N);
    for (int c = getc(fd); c!=EOF; c = getc(fd)) {
	if (isdigit(char(c))) {
	    ungetc(char(c), fd);
            int id;
            double vol, area;
            fscanf(fd, "%d %lf %lf", &id, &vol, &area);
	    printf("ID %d, vol %.2f, area %.2f\n", id, vol/1e6, area/1e4);
            if (name_map.find(id) == name_map.end() ) {
		printf("Size file lists ID %d not in name map!\n", id);
		return 42;
		}
            volumes[name_map[id]] = vol/1e6;  // 10nm per voxel
	    areas[name_map[id]]   = area/1e4;
	    }
	}
    }
if (fd != NULL)
    fclose(fd);
//-------------------------  Read the file of (potentially) incomplete bodies
vector<bool> incomplete(names.size(), true);  // start by assuming all are incomplete
bool csv = true;
//-------------------------  Read the file of neurotransmitter types
fd = fopen("nlist.txt", "r");
if (fd == NULL) {
    printf("No file 'nlist.txt' present\n");
    }
else {
    size_t nby;
    char *buffer = NULL;
    getline(&buffer, &nby, fd); // ignore the first line
    while (getline(&buffer, &nby, fd) != -1) {
        printf("Got '%s'\n", buffer); fflush(stdout);
	char *name = strtok(buffer, " ");
        char *tr = strtok(NULL, " ");	// get the second space separated token - the list of transmitters
        char *re = strtok(NULL, " ");   // list of receptors
        char *how_known = strtok(NULL, "\n");  // reason is the rest of the line
        char *how1 = strtok(how_known, ";");
        char *how2 = strtok(NULL,      "\n");
        NeuroTr t;
        t.name = strdup(name);
        t.trans = CommaList(tr);
        t.recep = CommaList(re);
        t.how_tr = how1 == NULL ? " " : strdup(how1);
        t.how_re = how2 == NULL ? " " : strdup(how2);
        NeuroTransmitters.push_back(t);
	}
    free(buffer);
    fclose(fd);
    }

//-------------------------  Read the synapse JSON file -------------------------------------------
// create a giant from-to matrix, with weights.
int M = names.size();
vector< vector<unsigned char> > from_to(M, vector<unsigned char>(M,0));
vector<int> mentioned(M,0);  // how often is each ID mentioned, either as a source or a target?
// Go through all the T-bars
//vector<TBar> tbs;
data = syn_root["data"];
for(unsigned int i=0; i<data.size(); i++) {
    //printf("getting entry %d\n", i);
    Json::Value v = data[i];
    // read the different things that can be in 'data'.  For now only {"T-bar","partner"} is known
    Json::Value t = v["T-bar"];
    Json::Value part = v["partners"];
    if (!t.isNull()) {
        Json::Value lo = t["location"];
        TBar tb;
	tb.status      = t["status"].asString();
	tb.confidence  = t["confidence"].asDouble();
	Json::Value co = t["convergent"];
        tb.convergent  = co.isNull() ? false : (co.asString() == "convergent");
        co             = t["flagged"];
        tb.flagged     = co.isNull() ? false : (co.asString() == "flagged");
	tb.body_id     = t["body ID"].asInt();
        co             = t["roi"];
        tb.roi         = co.isNull() ? 0 : (co.asString()[5] - '0');
        tb.gid         = 0;
	printf("TB convergent %d, flagged %d, roi %d\n", tb.convergent, tb.flagged, tb.roi);
        if (CT1s.find(tb.body_id) != CT1s.end())   // If any of the CT1s, replace by first
	    tb.body_id = *(CT1s.begin());
        if (name_map.find(tb.body_id) == name_map.end()) {
	    printf("No %d for TBar in map\n", tb.body_id);
	    continue;
	    }
	tb.pt.x = lo[0u].asDouble();
	tb.pt.y = lo[1u].asDouble();
        tb.z    = lo[2u].asInt();
	//printf("status is %s, confidence %f, loc %.1f,%lf,%d,  %d partners\n", 
	//tb.status.c_str(), tb.confidence, tb.pt.x, tb.pt.y, tb.z, part.size() );
	int from_id = name_map[tb.body_id];
	for(int j=0; j<part.size(); j++) {
	    Partner prt;
	    prt.confidence = part[j]["confidence"].asDouble();
	    prt.body_id    = part[j]["body ID"].asInt();
            Json::Value b  = part[j]["flagged"];
            prt.flagged    = b.isNull() ? false : b.asBool(); // if not present, make false
            printf("  Flag for psd is %d\n", prt.flagged);
	    if (CT1s.find(prt.body_id) != CT1s.end())   // If any of the CT1s, replace by first
		prt.body_id = *(CT1s.begin());
            if (name_map.find(prt.body_id) == name_map.end()) {
	        printf("No %d for partner in map\n", prt.body_id);
	        continue;
	        }
            int to_id = name_map[prt.body_id];
            if(from_to[from_id][to_id] == 255) {
		printf("Oh no!  Overflow! %d %d\n",tb.body_id, prt.body_id );
		}
	    if (KeepSynapse(tb.convergent)) {
                from_to[from_id][to_id]++;
                mentioned[from_id]++;
                mentioned[to_id]++;
	        }
	    Json::Value lo = part[j]["location"];
	    prt.pt.x = lo[0u].asDouble();
	    prt.pt.y = lo[1u].asDouble();
	    prt.z    = lo[2u].asInt();
            if (prt.body_id == tb.body_id) { // found an autapse
		char *name = names[name_map[prt.body_id]];
		if (!isdigit(name[0])) {
                    printf("Autapse:  %s Tbar at %d %d %d, psd at %d %d %d\n", name, 
		    int(tb.pt.x), int(tb.pt.y), tb.z, int(prt.pt.x), int(prt.pt.y), prt.z);
		    }
		}
	    //printf("   Partner %d, confidence %f. loc=%f %f %d\n", j, prt.confidence, prt.pt.x, prt.pt.y, prt.z);
	    tb.partners.push_back(prt);
	    }
	if (KeepSynapse(tb.convergent))
	    tbs.push_back(tb);
	}
    }
// If a cell is mentioned more than MENT_THRESHOLD times, but has no name, and it looks like a KC, make it a KC
const int MENT_THRESHOLD = 170;
int frag_id = 0;
for(int i=0; i<names.size(); i++) {
    if (mentioned[i] >= MENT_THRESHOLD && (isdigit(names[i][0]) || strncmp(names[i],"KC",2) == 0)) {
        // OK, look at total outs and ins, by partner count and total.
        int pre = 0, post = 0;       // sums
        int num_pre=0; int num_post=0;  // number of entries
        for(int y=0; y<names.size(); y++)
	    if (from_to[y][i] > 0) {
		pre += from_to[y][i];
                num_pre++;
	        }
        for(int x=0; x<names.size(); x++)
	    if (from_to[i][x] > 0) {
		post += from_to[i][x];
                num_post++;
	        }
        printf(" --- Name %s, pre %d/%d,   post %d/%d\n", names[i], pre, num_pre, post, num_post);
        if (isdigit(names[i][0]) &&
         70 <=  pre &&   pre <= 130 && 40 <=  num_pre &&   num_pre <= 80 &&
	 130 <= post && post <= 250 && 70 <= num_post && num_post <= 130) { //looks like a KC
	    char new_name[128];
            if (num_post >= 100)
	        sprintf(new_name, "KCs-guess-%d", frag_id++);
            else
	        sprintf(new_name, "KCc-guess-%d", frag_id++);
	    printf("Cell with name %s had %d mentions and is renamed to %s\n", names[i], mentioned[i], new_name);
	    names[i] = strdup(new_name);  // memory leak, but who cares
            }
	}
    }
// If a cell is mentioned more than MENT_THRESHOLD times, but has no name, make it "Unknown".
for(int i=0; i<names.size(); i++) {
    if (mentioned[i] >= MENT_THRESHOLD && isdigit(names[i][0])) {
	char new_name[128];
        sprintf(new_name, "Fragment-%d", frag_id++);
	printf("Cell with name %s had %d mentions and is renamed to %s\n", names[i], mentioned[i], new_name);
        names[i] = strdup(new_name);  // memory leak, but who cares
	}
    }
// If a cell is called "KC-a" look at its partners.  If all one of KC-c, KC-s, or KC-p, make it that type.
for(int changes=0; changes > 0; ){
    changes = 0;
    for(int i=0; i<names.size(); i++) {
	if (strcmp(names[i], "KC-any") != 0)
	    continue;
	int p=0;
	int c=0;
	int s=0;
	for(int k=0; k<names.size(); k++) {
	    if (from_to[i][k] == 0)
		continue;
	    if (strcmp(names[k], "KC-c") == 0) c++;
	    if (strcmp(names[k], "KC-s") == 0) s++;
	    if (strcmp(names[k], "KC-p") == 0) p++;
	    }
	for(int k=0; k<names.size(); k++) {
	    if (from_to[k][i] == 0)
		continue;
	    if (strcmp(names[k], "KC-c") == 0) c++;
	    if (strcmp(names[k], "KC-s") == 0) s++;
	    if (strcmp(names[k], "KC-p") == 0) p++;
	    }
	printf("GUESS c=%d s=%d p=%d\n",c, s, p);
        char *old = names[i];
	if (c > 10*(s+p)) {names[i] = strdup("KC-c"); changes++;}
	if (s > 10*(c+p)) {names[i] = strdup("KC-s"); changes++;}
	if (p > 10*(c+s)) {names[i] = strdup("KC-p"); changes++;}
        if (names[i] != old)
	    printf("Changed BodyId %d to %s\n", rev_map[i], names[i]);
	}
    printf("---GUESS: %d changes\n", changes);
    }
// Let's make a partitioning file for KCs
if (true) {
    FILE *fp = fopen("hm","w");             // One file for HMetis
    if (fp == NULL) {
	printf("Can't open hm\n");
	return 42;
	}
    vector<int>KCs;
    for(int i=0; i<names.size(); i++)
	if (strncmp(names[i],"KC",2) == 0)
	    KCs.push_back(i);
    printf("PART: There are %d Kenyon Cells\n", KCs.size() );
    // count the connections
    int count = 0;
    for(int i=0; i<KCs.size(); i++)
	for(int j=i+1; j<KCs.size(); j++)
	    count += (from_to[i][j] + from_to[j][i] > 0);
    // print the file.  With edge strengths
    fprintf(fp,"%d %d 1\n", count, KCs.size() );
    for(int i=0; i<KCs.size(); i++) 
	for(int j=i+1; j<KCs.size(); j++) {
	    int str = from_to[KCs[i]][KCs[j]] + from_to[KCs[j]][KCs[i]];
	    if(str > 0)
		fprintf(fp,"%d %d %d\n", str, i+1, j+1);
            }
    fclose(fp);
    fp = fopen("FixFile","w");
    if (fp == NULL) {
	printf("Can't open FixFile\n");
	return 42;
	}
    vector<int> part_count(2,0);
    for(int i=0; i<KCs.size(); i++) {
        int part = -1;
	if (part_count[0] < 10 && strcmp("KC-c", names[KCs[i]]) == 0) part = 0;
	if (part_count[1] < 10 && strcmp("KC-s", names[KCs[i]]) == 0) part = 1;
	fprintf(fp,"%d\n", part);
        if (part >= 0)
	     part_count[part]++;
	}
    fclose(fp);
    fflush(stdout);
    system("hmetis-2.0 hm 2 -ufactor=10 -fixed=\"FixFile\"");
    fflush(stdout);
    fp = fopen("hm.part.2","r");
    if (fp == NULL) {
	printf("Can't open hm.part.2\n");
	return 42;
        }
    for(int i=0; i<KCs.size(); i++) {
        int part;
	fscanf(fp, "%d", &part);
	//printf("Cell %s partition %d\n", names[KCs[i]], part);
	}
    fclose(fp);
    // Find the connected components
    vector<bool> used(KCs.size(), false);
    for(int i=0; i<KCs.size(); i++) {
	if (used[i])
	    continue;	// skip an already used
	printf("Component starting at %d\n", i);
	stack<int> look;
	look.push(i);
	while (!look.empty()) {
	    int j = look.top(); look.pop();
            if (used[j])
		continue;
	    used[j] = true;
	    for(int k=0; k<KCs.size(); k++) {
                 if (from_to[KCs[j]][KCs[k]] > 0) look.push(k);
                 if (from_to[KCs[k]][KCs[j]] > 0) look.push(k);
		 }
	    }
	 }
    // Do the same for Pajek
    // Is is better to uses 'arcs' (directed) or 'edges' (undirected)?
    fp = fopen("PajekTest.net","w");             // One file for HMetis
    if (fp == NULL) {
	printf("Can't open PajekTest.net\n");
	return 42;
	}
    fprintf(fp, "*vertices %d\n", KCs.size());
    for(int i=0; i<KCs.size(); i++)
	fprintf(fp, "%d \"%s-%d\"\n", i+1, names[KCs[i]], i+1);
    bool arcs = true;
    if (arcs) {
	fprintf(fp, "*arcs\n");
	for(int i=0; i<KCs.size(); i++)
	    for(int j=0; j<KCs.size(); j++) {
		int str = from_to[KCs[i]][KCs[j]];
		if(str > 0)
		    fprintf(fp,"%d %d %d\n", i+1, j+1, str);
		}
	}
    else { //edges
	fprintf(fp, "*edges\n");
	for(int i=0; i<KCs.size(); i++)
	    for(int j=i+1; j<KCs.size(); j++) {
		int str = from_to[KCs[i]][KCs[j]] + from_to[KCs[j]][KCs[i]];
		if(str > 0)
		    fprintf(fp,"%d %d %d\n", i+1, j+1, str);
		}
	}
    fclose(fp);
    vector<vector<Point3d> > rslts;
    rslts.push_back(ReadPajekFile("P-Eigen-2D.net"));
    rslts.push_back(ReadPajekFile("P-Eigen-3D.net"));
    rslts.push_back(ReadPajekFile("P-En-Fru-2D.net"));
    rslts.push_back(ReadPajekFile("P-En-Fru-3D.net"));
    rslts.push_back(ReadPajekFile("P-En-Kamada.net"));
    rslts.push_back(ReadPajekFile("P-MDS-2D.net"));
    rslts.push_back(ReadPajekFile("P-MDS-3D.net"));
    rslts.push_back(ReadPajekFile("P-VOS-2D.net"));
    rslts.push_back(ReadPajekFile("P-VOS-3D.net"));
    for(int i=0; i<KCs.size(); i++) {
        printf("%6d %8d %12s ", i, rev_map[KCs[i]], names[KCs[i]]);
        vector<char> cs;
	for(int k=0; k<rslts.size(); k++) {
	    char c = FindMembership(i, rslts[k], KCs, names);
	    printf("%c", c);
            cs.push_back(c);
	    }
	printf("\n");
        char *name = names[KCs[i]];
        if (strcmp(name,"KC-s") == 0 || strcmp(name,"KC-p") == 0 || strcmp(name, "KC-c") == 0) {
	    // already classified.  See if consistent.  In best case, get all upper case of same type
	    // (upper case means 2 closest both agree).
	    char cl = toupper(name[3]);
            int agree = 0, disagree = 0;
            for(int k=0; k<cs.size(); k++) {
		agree    += (cs[k] == cl);
		disagree += (cs[k] != cl);
		}
	    if (disagree > 1) {
		printf("--------Susicious classification.  Body ID %d is labeled %s but best matches are ", rev_map[KCs[i]], name);
                for(int k=0; k<cs.size(); k++)
		    printf("%c",cs[k]);
		printf("\n");
		}
	    }
	else if (strcmp(name,"KC-any") == 0) {
	    char cl = tolower(cs[0]);
	    int agree = 0, disagree = 0;
	    for(int k=0; k<cs.size(); k++) {
		agree    += (tolower(cs[k]) == cl);
		disagree += (tolower(cs[k]) != cl);
		}
	    //printf("agree %d, disagree %d, size %d\n", agree, disagree, cs.size() );
	    if (agree == cs.size()) {
		printf("--------It's unanimous.  Body ID %d is labeled %s but looks like KC-%c\n", rev_map[KCs[i]], name, cl);
		}
	    }
        }
    // make another array of types
    FindCommonInputs("MBON-19", KCs, names, from_to, rev_map);
    FindCommonInputs("MBON-14", KCs, names, from_to, rev_map);
    FindCommonInputs("MBON-07", KCs, names, from_to, rev_map);
    }
   
// Look for partners for convergent synpases.
const int TOO_FAR = 500;  // surely not more than 500 pixels away
vector<vector<int> > groups;
int group_id = 1;
for(int i=0; i<tbs.size(); i++) {
    if (!tbs[i].convergent)
	continue;
    if (tbs[i].gid > 0)
	continue;
    // OK, it's convergent and not yet part of a group.
    const double CLOSE = 50;  // 30 pixels, for now
    stack<int> st;
    tbs[i].gid = group_id;
    st.push(i);
    printf("\nGroup %d, adding %d\n", group_id, i);
    vector<int> gp;
    while (!st.empty()) {
	int curr = st.top();
	st.pop();
	printf(" Just popped %d\n", curr);
        gp.push_back(curr);
	double best = 1e5;
        int old = st.size();
        for(int j=0; j<tbs[curr].partners.size(); j++) {
	    for(int k=0; k<tbs.size(); k++) {
		if (tbs[k].gid > 0)
		    continue;
		if (abs(tbs[curr].z - tbs[k].z) > TOO_FAR) continue;
		for(int m=0; m<tbs[k].partners.size(); m++) {
		    double d = PartnerDist(tbs[curr].partners[j], tbs[k].partners[m]);
                    best = min(best, d);
		    if (d <= CLOSE) {
			printf("Tbar %d, [%d] <- %.2f -> Tbar %d, [%d] \n", i, j, d, k, m);
			tbs[k].gid = group_id;
			st.push(k);
			printf(" Pushing %d\n", k);
		        break;
			}
		    }
		}
	    }
        if (st.size() == old)
	    printf("None found within range.  Best = %.2f\n", best);
	}
    groups.push_back(gp);
    group_id++;
    }
// Print them out.
for(int i=0; i<groups.size(); i++) {
    printf("\nGroup %d\n", i);
    map<int,int> m;
    for(int j=0; j<groups[i].size(); j++) {
        int t = groups[i][j];
	printf("  Tbar id %8d%c%c -> ", tbs[t].body_id, tbs[t].convergent ? '*': ' ', tbs[t].flagged ? '+' : ' ');
	for(int k=0; k < tbs[t].partners.size(); k++) {
            int n = tbs[t].partners[k].body_id;
	    printf("%8d%c ", n, tbs[t].partners[k].flagged ? '+' : ' ');
            if (m.find(n) == m.end() )
		m[n] = 1;
	    else
		m[n] = m[n]+1;
	    }
	printf("\t(%5d, %5d, %5d)", int(tbs[t].pt.x), int(tbs[t].pt.y), tbs[t].z);
	printf("\n");
	}
    map<int,int>::iterator mi;
    for(mi = m.begin(); mi != m.end(); mi++)
	printf("%9d %2d%c\n", mi->first, mi->second, mi->second == groups[i].size() ? '*' : ' ');
   }
// Zipursky code here ------------------------------------------------------------
int tiL2 = TypeIndex("L2");
printf("Zipurksy, looking for %d\n", tiL2);
for(unsigned int i=0; i<data.size(); i++) {
    //printf("getting entry %d\n", i);
    Json::Value v = data[i];
    // read the different things that can be in 'data'.  For now only {"T-bar","partner"} is known
    Json::Value t = v["T-bar"];
    Json::Value part = v["partners"];
    if (!t.isNull()) {
        Json::Value lo = t["location"];
        TBar tb;
	tb.status      = t["status"].asString();
	tb.confidence  = t["confidence"].asDouble();
	tb.body_id     = t["body ID"].asInt();
        if (name_map.find(tb.body_id) == name_map.end()) {
	    printf("No %d for TBar in map\n", tb.body_id);
	    continue;
	    }
	char *name = names[name_map[tb.body_id]];
        if (TypeIndex(name) != tiL2)
	    continue;
	tb.pt.x = lo[0u].asDouble();
	tb.pt.y = lo[1u].asDouble();
        tb.z    = lo[2u].asInt();
	printf(" Tb-bar on %s, partners [", name );
        vector<string> pnames;
	for(int j=0; j<part.size(); j++) {
	    Partner prt;
	    prt.confidence = part[j]["confidence"].asDouble();
	    prt.body_id    = part[j]["body ID"].asInt();
            if (name_map.find(prt.body_id) == name_map.end()) {
	        printf("No %d for partner in map\n", prt.body_id);
	        continue;
	        }
            name = names[name_map[prt.body_id]];
            if (isdigit(name[0]))
		continue;
	    Json::Value lo = part[j]["location"];
	    prt.pt.x = lo[0u].asDouble();
	    prt.pt.y = lo[1u].asDouble();
	    prt.z    = lo[2u].asInt();
	    pnames.push_back(string(name));
	    }
	sort(pnames.begin(), pnames.end());
        int ndup = 0;
	for(int j=0; j<pnames.size(); j++) {
	    printf("%12s", pnames[j].c_str());
            if (j > 0 && pnames[j].compare(pnames[j-1]) == 0)
		ndup++;
	    }
	printf(" ] %d dup\n", ndup);
	}
    }
printf("End Zipursky\n");
// End Zipursky ------------------------------------------------------------------
// Percentage complete code
printf("Completeness by PSDs identified\n");
for(int from=0; from<names.size(); from++) {
    if (isdigit(names[from][0]))
	continue;
    int known = 0; int unknown = 0;
    for(int to=0; to<names.size(); to++) {
	int str = from_to[from][to];
	if (str == 0) 
	    continue;
	if (isdigit(names[to][0]))
	    unknown += str;
	else
	    known += str;
	}
    printf(" Cell %15s: %4d known, %4d unknown, %.2f%% complete\n", names[from], known, unknown, known/double(known+unknown)*100.0);
    }
printf("End of Completeness\n");
// End of completeness -----------------------------------------------------------
bool slow = true;
if (slow) {
printf("Unknown cells:\n");
vector<sort_down> sd;
for(int i=0; i<names.size(); i++) {
    if (strncmp(names[i], "Frag", 4) != 0)  // only try the fragments
	continue;
    // find output, known and unknown, and inputs, known and unknown
    int ink = 0; int inu = 0, outk = 0, outu = 0;
    for(int to=0; to<names.size(); to++) {
	int str = from_to[i][to];
	if (str == 0) 
	    continue;
	if (isdigit(names[to][0]))
	    outu += str;
	else
	    outk += str;
	}
    for(int from=0; from<names.size(); from++) {
	int str = from_to[from][i];
	if (str == 0) 
	    continue;
	if (isdigit(names[from][0]))
	    inu += str;
	else
	    ink += str;
	}
    printf(" Cell %15s: %4d known in, %4d unknown in, %4d known out, %4d unknown out, %4d total\n", names[i], ink, inu, outk, outu, ink+inu+outk+outu);
    fflush(stdout);
    sort_down s(i, ink + outk);
    sd.push_back(s);
    }
printf("End of Unknown.  %d candidates\n", sd.size() );
fflush(stdout);
sort(sd.begin(), sd.end());
for(int j=0; j<25 && j < sd.size(); j++) {
    int i = sd[j].index;
    printf("Cell %s has %d known connections\n", names[i], int(sd[j].val) );
    vector<int> ins_by_type(NumTypes,0);
    vector<int> out_by_type(NumTypes,0);
    vector<int> str_ins_by_type(NumTypes, 0);
    vector<int> str_out_by_type(NumTypes, 0);
    for(int to=0; to<names.size(); to++) {
	int str = from_to[i][to];
	if (str == 0) 
	    continue;
	if (!isdigit(names[to][0])) {
	    printf("   To %15s, strength %d\n", names[to], str);
            int ty = TypeIndex(names[to]);
            out_by_type[ty]++;
            str_out_by_type[ty] += str;
	    }
	}
    for(int from=0; from<names.size(); from++) {
	int str = from_to[from][i];
	if (str == 0) 
	    continue;
	if (!isdigit(names[from][0])) {
	    printf(" From %15s, strength %d\n", names[from], str);
            int ty = TypeIndex(names[from]);
            ins_by_type[ty]++;
            str_ins_by_type[ty] += str;
	    }
	}
    for(int k=0; k<NumTypes; k++) {
        if (ins_by_type[k] > 0 || out_by_type[k] > 0) {
	    printf("  Type %15s, ins %d cells, total str %d,  outs %d cells, total str %d\n", CellTypes[k].name,
	     ins_by_type[k], str_ins_by_type[k], out_by_type[k], str_out_by_type[k] );
	    }
	}
    }
}
// End of completeness -----------------------------------------------------------
printf("Read %d Tbars\n", tbs.size());
// Plot maps of Synapse spots
// These should look up the synapse IDS from the names, but we cheat here.
//MapXYs(tbs, 568121, "bf", name_map, names);
//MapXYs(tbs, 571372, "du", name_map, names);
//MapXYs(tbs, 597861, "fb", name_map, names);
//MapXYs(tbs, 598856, "ud", name_map, names);
vector<double> Mi1x;   // keep track of all 'Mi1 H' specifically, to fit column coordinates.
vector<double> Mi1y;
vector<double> Mi1z;
// OK, create the set to average the location of each chunk.   This is the only use of the tsb vector, so this could
// be done in the above loop.  Also, since we now create 'ss' we could get rid of the huge matrix.
// Also, create a matrix of which cell types have synapses in each layer.
// Also, create a T-bar count for each cell.
vector<vector<int> > s_vs_ml(NumTypes, vector<int>(11, 0));  // synapses in each layer, by type
vector<int> Tbar_count(names.size(), 0);
vector<vector<int> > Tbar_count_by_layer(names.size(), vector<int>(11,0));
set<where> Tbar_by_pair;   // a 'where' data structure is overkill for this, but should work and use
			   // much less memory than a matrix
for(int i=0; i<tbs.size(); i++) {
    Tbar_count[name_map[tbs[i].body_id]]++;
    int med_lay = XYZtoLayer(int(tbs[i].pt.x), int(tbs[i].pt.y), tbs[i].z);
    Tbar_count_by_layer[name_map[tbs[i].body_id]][0]++;        // 0 never occurs normally, so use it for sum of all layers
    Tbar_count_by_layer[name_map[tbs[i].body_id]][med_lay]++;
    for(int j=0; j<tbs[i].partners.size(); j++) {
        int id = name_map[tbs[i].partners[j].body_id];
        where local(name_map[tbs[i].body_id], id, med_lay);
        //printf("## Insert %d %d %d\n", local.from_id, local.to_id, med_lay);  // convert to microns
	pair<set<where>::iterator, bool>  ret = ss.insert(local);
        ret.first->sx += int(tbs[i].pt.x);
        ret.first->sy += int(tbs[i].pt.y);
        ret.first->sz += int(tbs[i].z);
        ret.first->n++;
        int ti = TypeIndex(names[id]);
        s_vs_ml[ti][med_lay]++;
	}
    if (name_map[tbs[i].body_id] == Mi1) {
	Mi1x.push_back(tbs[i].pt.x);
	Mi1y.push_back(tbs[i].pt.y);
	Mi1z.push_back(tbs[i].z);
	}
    // create a unique list of partners
    vector<int>partners;
    for(int j=0; j<tbs[i].partners.size(); j++)
	partners.push_back(tbs[i].partners[j].body_id);
    sort(partners.begin(), partners.end());
    for(int j=0; j<partners.size(); j++) {
	if (j == 0 || partners[j] != partners[j-1]) { // it's unique
	    int id = name_map[partners[j]];
	    where local(name_map[tbs[i].body_id], id, 0);  // use medulla layer 0
	    //printf("## Insert %d %d %d\n", local.from_id, local.to_id, med_lay);  // convert to microns
	    pair<set<where>::iterator, bool>  ret = Tbar_by_pair.insert(local);
	    ret.first->n++;
	    }
	}
    }
// ----------------- For named pairs, compute strength/Tbar ratio
printf("Strength/TBar ratio\n");
for(int from=0; from<names.size(); from++) {
    if (isdigit(names[from][0]))
	continue;
    for(int to=0; to<names.size(); to++) {
	int str = from_to[from][to];
	if (str == 0 || isdigit(names[to][0]))
	    continue;
        where local(from, to, 0);
        set<where>::iterator look = Tbar_by_pair.find(local);
        if (look == Tbar_by_pair.end()) // should always find at least one.
	    printf("What?? can't find %s %s %d %d\n", names[from], names[to], from, to);
	else
	    printf("  From %11s, to %11s, ratio=%.3f (%4d conn, %4d Tbars)\n", names[from], names[to], double(str)/look->n, str, look->n);
	}
    }
printf("End of ratio\n");
map<string, double> from_to_ratio;
map<string, double> from_to_weight;
// Now, let's do this again, assuming the unknown partners are picked at random from the same cell's known connections in the same layer.
set<where> str_by_pair;   // a 'where' data structure is overkill for this, but should work and use
			  // little memory
{
// First index is the ID into names[], second is the medulla layer, third is the list of body_ids;
vector<vector<vector<int> > > ctab(names.size(), vector<vector<int> >(11));
for(int i=0; i<tbs.size(); i++) {
    int T_id = name_map[tbs[i].body_id];
    if (isdigit(names[T_id][0]))
	continue;
    int med_lay = XYZtoLayer(int(tbs[i].pt.x), int(tbs[i].pt.y), tbs[i].z);
    for(int j=0; j<tbs[i].partners.size(); j++) {
        int id = tbs[i].partners[j].body_id;
        if (isdigit(names[name_map[id]][0]))
	    continue;
        ctab[T_id][med_lay].push_back(id);
	}
    }
// Now, for all un-named partners, set according to layer specific table.
for(int i=0; i<tbs.size(); i++) {
    int T_id = name_map[tbs[i].body_id];
    if (isdigit(names[T_id][0]))
	continue;
    int med_lay = XYZtoLayer(int(tbs[i].pt.x), int(tbs[i].pt.y), tbs[i].z);
    for(int j=0; j<tbs[i].partners.size(); j++) {
        int id = name_map[tbs[i].partners[j].body_id];
        if (isdigit(names[id][0])) {
            int N = ctab[T_id][med_lay].size();
	    if (N == 0) {
		printf("No named connections for %s, layer M%d\n", names[T_id], med_lay);
	        tbs[i].partners[j].fake_id = tbs[i].partners[j].body_id;   // so it still could un-named, worst case
	        }
            else
	        tbs[i].partners[j].fake_id = ctab[T_id][med_lay][random() % N];
            }
	else  // it has a real name.  Just copy it
	    tbs[i].partners[j].fake_id = tbs[i].partners[j].body_id;
	}
    }
// Now, redo the calculation using the fake_ids.  All the strengths will change too, so redo them
Tbar_by_pair.clear();
for(int i=0; i<tbs.size(); i++) {
    int T_id = name_map[tbs[i].body_id];
    if (isdigit(names[T_id][0]))
	continue;
    // Find the new strengths.  These do not care about uniqueness
    for(int j=0; j<tbs[i].partners.size(); j++) {
        int id = tbs[i].partners[j].fake_id;
        if (isdigit(names[name_map[id]][0]))
	    continue;
	where local(T_id, name_map[id], 0);
	pair<set<where>::iterator, bool>  ret = str_by_pair.insert(local);
	ret.first->n++;
	}
    // Now T-bars.  these *do* care about uniqueness.
    vector<int>partners;
    for(int j=0; j<tbs[i].partners.size(); j++)
	partners.push_back(tbs[i].partners[j].fake_id);
    sort(partners.begin(), partners.end());
    for(int j=0; j<partners.size(); j++) {
	if (j == 0 || partners[j] != partners[j-1]) { // it's unique
	    int id = name_map[partners[j]];
	    where local(name_map[tbs[i].body_id], id, 0);  // use medulla layer 0
	    //printf("## Insert %d %d %d\n", local.from_id, local.to_id, med_lay);  // convert to microns
	    pair<set<where>::iterator, bool>  ret = Tbar_by_pair.insert(local);
	    ret.first->n++;
	    }
	}
    }
printf("2:Strength/TBar ratio\n");
for(int from=0; from<names.size(); from++) {
    if (isdigit(names[from][0]))
	continue;
    for(int to=0; to<names.size(); to++) {
	if (isdigit(names[to][0]))
	    continue;
        where local(from, to, 0);
	int str = 0;
        set<where>::iterator look = str_by_pair.find(local);
        if (look != str_by_pair.end())
	    str = look->n;
        if (str == 0)
	    continue;
        look = Tbar_by_pair.find(local);
        if (look == Tbar_by_pair.end()) // should always find at least one.
	    printf("What?? can't find %s %s %d %d\n", names[from], names[to], from, to);
	else {
	    printf("  From %11s, to %11s, ratio=%.3f (%4d conn, %4d Tbars)\n", names[from], names[to], double(str)/look->n, str, look->n);
            string name = string(BaseName(names[from])) + ":" + string(BaseName(names[to]));
	    printf("Inserting %s\n", name.c_str() );
            if (from_to_ratio.find(name) != from_to_ratio.end()) {
	        from_to_ratio[name]  = from_to_ratio[name]  + double(str);
                from_to_weight[name] = from_to_weight[name] + look->n;
		}
	    else {
	        from_to_ratio[name]  = double(str);
                from_to_weight[name] = look->n;
		}
	    }
	}
    }
printf("2:End of ratio\n");
}
printf("There are %d unique ID pairs\n", ss.size());
// Print a table of cell types and synapse counts
printf("Cell type     M1    M2    M3    M4    M5    M6    M7    M8    M9   M10\n");
printf("----------  ----  ----  ----  ----  ----  ----  ----  ----  ----  ----\n");
vector<int> HowManyTypes(11,0);
for(int i=0; i<NumTypes; i++) {
    printf("%10s", CellTypes[i].name);
    for(int j=1; j<=10; j++) {
	printf("%6d", s_vs_ml[i][j] );
        if (s_vs_ml[i][j] >= 10)
	    HowManyTypes[j]++;
	}
    printf("\n");
    }
printf("----------  ----  ----  ----  ----  ----  ----  ----  ----  ----  ----\n");
printf("Types >=10");
for(int j=1; j<=10; j++)
    printf("%6d", HowManyTypes[j]);
printf("\n");
// Print the names of all cells with Tbars
printf("\nTbar counts\n");
for(int i=0; i<names.size(); i++) {
    if (Tbar_count[i] > 0)
	printf("Cell %12s has %d Tbars\n", names[i], Tbar_count[i]);
    }
// Write a file of overall type information.   This is a set of 10 matrices, each NumTypes x Numtypes.
// One has the total connection count, then one for each layer.
{
vector<vector<vector<int> > > from_to_by_type(NumTypes, vector<vector<int> > (NumTypes, vector<int>(11, 0)));
vector<int> TBs(NumTypes, 0);  // count TBars
for(int i=0; i<tbs.size(); i++) {
    int T_id = name_map[tbs[i].body_id];
    if (isdigit(names[T_id][0]))
	continue;
    int med_lay = XYZtoLayer(int(tbs[i].pt.x), int(tbs[i].pt.y), tbs[i].z);
    int T_type = TypeIndex(names[T_id]);
    TBs[T_type]++;
    for(int j=0; j<tbs[i].partners.size(); j++) {
        int id = name_map[tbs[i].partners[j].body_id];
        if (isdigit(names[id][0]))
	    continue;
        int P_type = TypeIndex(names[id]);
        from_to_by_type[T_type][P_type][med_lay]++;  // by medulla layer
        from_to_by_type[T_type][P_type][0]++;        // layer 0, not otherwise used, gets the sum
	}
    }
// Now print the matrices in comma separated form.
for(int k=0; k<11; k++) {
    char name[256];
    sprintf(name, "ConnByType.%d.csv", k);
    FILE *f = fopen(name, "w");
    if (f == NULL) {
	printf("Cound not open %s\n", name);
	return 42;
	}
    // write the top line
    for(int i=0; i<NumTypes; i++)
	fprintf(f, ",%s", CellTypes[i].name);
    fprintf(f, "\n");
    for(int from = 0; from < NumTypes; from++) {
        // write each line
	fprintf(f, "%s", CellTypes[from].name);
	for(int to=0; to < NumTypes; to++)
	    fprintf(f, ",%d", from_to_by_type[from][to][k] );
	fprintf(f, "\n");
	}
    fclose(f);
    }
// print TBars for Michael
printf("TBars for Michael\n");
for(int i=0; i<NumTypes; i++)
    printf("%15s %4d\n", CellTypes[i].name, TBs[i] );
printf("End Michael\n");
}
// Write a file for clustering attempts.  Compute a feature vector for each named cell.
printf("Clustering test\n");
fflush(stdout);
for(int i=0; i<names.size(); i++) {
    if (incomplete[i] || isdigit(names[i][0]))
	continue;
    vector<int>  ins(11,0);  // need indices 1-10
    vector<int> outs(11,0);
    for(int j=0; j<names.size(); j++) {
	for(int k=1; k<=10; k++) {
	    where look(i, j, k);
            set<where>::iterator it = ss.find(look);
            if (it != ss.end())
                outs[k] += it->n;
	    where look2(j,i,k);
            it = ss.find(look2);
            if (it != ss.end())
                ins[k] += it->n;
	    }
	}
    printf("%12s :", names[i]);
    fflush(stdout);
    for(int k=1; k<=10; k++)
	printf("%5d", ins[k]);
    printf("    ");
    for(int k=1; k<=10; k++)
	printf("%5d", outs[k]);
    printf("%8.2f %8.2f\n", volumes[i], areas[i]);
    fflush(stdout);
    }

// Calculate tilt of the columns
if (Mi1x.size() > 1) {
    printf("Estimating from %d synapses on 'Mi1 H':\n", Mi1x.size() );
    Fitab xfitter(Mi1z,Mi1x);
    printf("  X as a function of Z %f %f\n", xfitter.a, xfitter.b);
    Fitab yfitter(Mi1z,Mi1y);
    printf("  Y as a function of Z %f %f\n", yfitter.a, yfitter.b);
    }
else {
    printf("No Mi1 H??\n");
    //return 42;
    }

// Look for all the overlaps
fd = fopen("gimme_overlaps","w");  // write a file of scripts to get volumes and areas
if (fd == NULL) {
    printf("Could not open script file 'gimme_overlap'\n");
    return 42;
    }
fprintf(fd, "#!/bin/sh\n");
fprintf(fd, "rm overlaps\n");
for(int i=0; i<names.size(); i++) {
    if (isdigit(names[i][0]))
	continue;
    fprintf(fd, "curl -X GET http://emdata2:9000/api/node/fa9/graph/neighbors/%d >tmp.json\n", rev_map[i]);
    fprintf(fd, "python summary.py >>overlaps\n");
    }
fclose(fd);
system("chmod +x gimme_overlaps");

// Now read the overlaps, presumably from a previous run
set<overlap> Overlaps;
fd = fopen("overlaps","r");  // write a file of scripts to get volumes and areas
if (fd == NULL) {
    printf("Could not open file 'overlaps' for read\n");
    }
else {
    int i = 0;
    for(; ; i++) {
	int i1, i2;
        double d;
        if (fscanf(fd, "%d %d %lf", &i1, &i2, &d) != 3)
	    break;
        d = d/10000.0;   /// convert to square microns
	if (i2 < i1)
	    {int swap = i1; i1 = i2; i2 = swap;}
        overlap o(i1, i2, d);
	Overlaps.insert(o);
	}
    printf("Read %d overlaps\n", i);
    vector<sort_overlap> sort_over(Overlaps.size());
    int j=0;
    for(set<overlap>::iterator i = Overlaps.begin(); i != Overlaps.end(); i++)
	sort_over[j++].it = i;
    sort(sort_over.begin(), sort_over.end() );
    printf("Sorted!\n");
    for(int i=0; i<sort_over.size(); i++) {
        int id1 = name_map[sort_over[i].it->id1];
        int id2 = name_map[sort_over[i].it->id2];
        char *name1 = names[id1];
        char *name2 = names[id2];
        int n1 = from_to[id1][id2];
        int n2 = from_to[id2][id1];
        if (strcmp(name1,name2) > 0)
	    {char *tmp = name1; name1 = name2; name2 = tmp; int t = n1; n1 = n2; n2 = t;}
        int sum = from_to[id1][id2] + from_to[id2][id1];
        double area = sort_over[i].it->area;
	 printf("%9s %9s %10.4f - %3d %3d %3d  - %5.2f %5.2f %5.2f per square micron\n", name1, name2, sort_over[i].it->area, 
	 n1, n2, sum, n1/area, n2/area, sum/area);
	}
    return 42;
    }

// Calculate an graph for simulation.
vector<bool> wanted(NumTypes, false);
vector<bool> trace (NumTypes, true);
wanted[TypeIndex("Tm3")] = true;
wanted[TypeIndex("Mi1")] = true;
wanted[TypeIndex("Mi4")] = true;
wanted[TypeIndex("Mi9")] = true;
wanted[TypeIndex( "L1")] = true; trace[TypeIndex("L1")] = false;
wanted[TypeIndex( "L2")] = true; trace[TypeIndex("L2")] = false;
wanted[TypeIndex( "L3")] = true; trace[TypeIndex("L3")] = false;
wanted[TypeIndex( "L4")] = true; trace[TypeIndex("L4")] = false;
wanted[TypeIndex( "L5")] = true; trace[TypeIndex("L5")] = false;
wanted[TypeIndex( "R7")] = true; trace[TypeIndex("R7")] = false;
wanted[TypeIndex( "R8")] = true; trace[TypeIndex("R8")] = false;
wanted[TypeIndex( "C2")] = true; trace[TypeIndex("C2")] = false;
const int threshold = 3;

int root = 568121;
CircuitBacktrace(name_map[root], from_to, names, wanted, trace, threshold);
root = 597861;
CircuitBacktrace(name_map[root], from_to, names, wanted, trace, threshold);

// sort the names by their assigned color.  (Better than names, which sort oddly in Mi1, Mi15, Mi1b, etc.)
// Since columns a, c,e, and f are pale, and others yellow ometidia , we want them to sort together.   So map, then unmap.
vector<sortname> sn(M);
for(int k=0; k<names.size(); k++) {
    //              a->  b->  c->  d->  e->  f->
    char lmap[6] = {'C', 'A', 'D', 'B', 'E', 'F'};
    sn[k].index = k;
    // Normally, we sort by color.  This is an integer 0..2^24-1.  All cells of the same type have the same color,
    // so this sorts by cell type.  Next we sort by name - this puts the columns A-F at the top of the table.
    // But for non-columnar cells, columns make no sense, so we'd like to sort by size.  We can do this by
    // adding a small (less than 1.0) fraction to the color.  Since color overrides name, this will sort them
    // correctly.  And for columnar cells, color will still be identical, so sort will go to second key, name.
    sn[k].color = color(names[k]);
    if ((!IsColumnar(names[k])) && volumes[k] > 0 && volumes[k] < 1e6) {
        printf("Tweaking %s by %f\n", names[k], volumes[k]/1e6);
	sn[k].color += 1-volumes[k]/1e6;
        }

    sn[k].name  = names[k];
    char *last = names[k] + strlen(names[k]) -1;
    if (*last >= 'A' && *last <= 'F')
	*last = lmap[*last - 'A'];
    }
sort(sn.begin(), sn.end() );

// reverse the mapping
for(int k=0; k<names.size(); k++) {
    //              a->  b->  c->  d->  e->  f->
    char lmap[6] = {'B', 'D', 'A', 'C', 'E', 'F'};
    char *last = names[k] + strlen(names[k]) -1;
    if (*last >= 'A' && *last <= 'F')
	*last = lmap[*last - 'A'];
    }

vector<sort_down> sa(M);  // create a sort array
int cc = 0;            // current color

// While we are printing, we will collect the data for the matrix equivalent to the one in the Single Column paper
//
const int NS = 27;  // Number of type in single column paper
vector<vector<double> > scmat(NS, vector<double>(NS, 0.0));    // all cells in reconstruction
vector<vector<double> > cen7mat(NS, vector<double>(NS, 0.0));  // Just in the central 7 columns
vector<vector<double> > homemat(NS, vector<double>(NS, 0.0));  // Just the home column; should be closest to previous paper
vector<vector<double> >   hhmat(NS, vector<double>(NS, 0.0));  // Just the home-home connections
vector<vector<double> > palemat(NS, vector<double>(NS, 0.0));  // for pale-pale
vector<vector<double> > yellmat(NS, vector<double>(NS, 0.0));  // yellow-yellow
vector<vector<double> > sc_palemat(NS, vector<double>(NS, 0.0)); // Within column, in a pale column
vector<vector<double> > sc_pale_sd(NS, vector<double>(NS, 0.0)); // Within column, in a pale column
vector<vector<double> > sc_yellmat(NS, vector<double>(NS, 0.0)); // Within column, in a yellow column
vector<vector<double> > sc_yell_sd(NS, vector<double>(NS, 0.0)); // Within column, in a yellow column

vector<vector<int > > out_rows;   // the actual rows
vector<vector<int > >  in_rows;   // the actual rows
vector<vector<int > >   out_rows_dat;
vector<vector<int > >    in_rows_dat;
vector<vector<int > >   out_rows_recip;
vector<vector<int > >    in_rows_recip;
vector<int>  in_row_sums;          // total numnber of inputs in each row
vector<int> out_row_sums;
vector<double>  in_row_area;
vector<double>  in_row_vols;
vector<double> out_row_area;
vector<double> out_row_vols;
vector<bool>  in_row_inc;	   // does each row represent an incomplete cell?
vector<bool> out_row_inc;
vector<int> in_row_ids;       // ids of the rows
vector<int>out_row_ids;       // ids of the rows

// For each cell type, create a vector of all the types connected to.
vector<MeanStd> out_stats(NumTypes, MeanStd());  // stats for all connections (not just ones in table)
vector<MeanStd> out_7col (NumTypes, MeanStd());  // stats for all connections in the main 7 columns
vector<MeanStd> out_pale (NumTypes, MeanStd());  // stats for all connections in the pale columns
vector<MeanStd> out_yell (NumTypes, MeanStd());  // stats for all connections in the yellow columns
vector<MeanStd> out_home (NumTypes, MeanStd());  // stats for just the home column
vector<MeanStd> out_hh   (NumTypes, MeanStd());  // stats for just home-home connections
vector<MeanStd>  in_stats(NumTypes, MeanStd());  // stats for all connections (not just ones in table)
vector<MeanStd>  in_7col (NumTypes, MeanStd());  // stats for all connections in the main 7 columns
vector<MeanStd>  in_pale (NumTypes, MeanStd());  // stats for all connections (pale,pale)
vector<MeanStd>  in_yell (NumTypes, MeanStd());  // stats for all connections (yellow, yellow) columns
vector<MeanStd>  out_sc_pale (NumTypes, MeanStd());  // stats for within column (ie. b-b) in pale columns
vector<MeanStd>   in_sc_pale (NumTypes, MeanStd());  // stats for within column (ie. b-b) in pale columns
vector<MeanStd>  out_sc_yell (NumTypes, MeanStd());  // stats for within column (ie. b-b) in pale columns
vector<MeanStd>   in_sc_yell (NumTypes, MeanStd());  // stats for within column (ie. b-b) in pale columns

vector<vector<int> >  in_cross_col(8, vector<int>(8,0));  // from one column to another
vector<vector<int> > out_cross_col(8, vector<int>(8,0));
const int NUM_SHOWN = 1000;
for(int ii=0; ii<M+1; ii++) {  // One extra time, to print last table
    int i;
    if (ii != M) {
        i = sn[ii].index;
        if (isdigit(names[i][0]))
	    continue;
        }

    // new color?  then print the old table, if any, then start new table
    if (ii == M || color(names[i]) != cc) {
        if (in_row_ids.size()> 0) 
	    write_html(names, rev_map, in_row_ids,  "in",  in_rows,   in_rows_dat,  in_rows_recip, in_row_sums, in_row_vols, in_row_area, in_row_inc,
            in_stats, in_7col, in_pale, in_yell, in_cross_col, Tbar_count);
        if (out_row_ids.size() > 0) {
	    write_html(names, rev_map, out_row_ids, "out", out_rows, out_rows_dat, out_rows_recip, out_row_sums, out_row_vols, out_row_area, out_row_inc,
            out_stats, out_7col, out_pale, out_yell, out_cross_col, Tbar_count);
            // add the entries of the current table to the single column (from,to) matrix
            int ti = TypeIndex(names[out_row_ids[0]]);
            if (CellTypes[ti].sc_index >= 0) { // it was in the SC paper
                int row = CellTypes[ti].sc_index;
                for(int k=0; k<NumTypes; k++) {
                    if (row == 0)
			printf("First row.  k=%d->%d.  %d entries, %f mean\n", k, CellTypes[k].sc_index, out_stats[k].HowMany(), out_stats[k].Mean() );
	            if(out_stats[k].HowMany() > 0 && CellTypes[k].sc_index >= 0)
		        scmat[row][CellTypes[k].sc_index] = out_stats[k].Mean()*double(out_stats[k].HowMany())/out_row_ids.size();
	            if(out_7col[k].HowMany() > 0 && CellTypes[k].sc_index >= 0)
		        cen7mat[row][CellTypes[k].sc_index] = out_7col[k].Mean()*double(out_7col[k].HowMany())/7.0;
	            if(out_pale[k].HowMany() > 0 && CellTypes[k].sc_index >= 0)
		        palemat[row][CellTypes[k].sc_index] = out_pale[k].Mean()*double(out_pale[k].HowMany())/NPALE;
	            if(out_yell[k].HowMany() > 0 && CellTypes[k].sc_index >= 0)
		        yellmat[row][CellTypes[k].sc_index] = out_yell[k].Mean()*double(out_yell[k].HowMany())/NYELLOW;
	            if(out_home[k].HowMany() > 0 && CellTypes[k].sc_index >= 0)
		        homemat[row][CellTypes[k].sc_index] = out_home[k].Mean()*double(out_home[k].HowMany());
	            if(out_hh  [k].HowMany() > 0 && CellTypes[k].sc_index >= 0)
		        hhmat  [row][CellTypes[k].sc_index] = out_hh[k].Mean()*double(out_hh[k].HowMany());
	            if(out_sc_pale[k].HowMany() > 0 && CellTypes[k].sc_index >= 0) {
		        sc_palemat[row][CellTypes[k].sc_index] = out_sc_pale[k].Mean()*double(out_sc_pale[k].HowMany())/NPALE;
                        sc_pale_sd[row][CellTypes[k].sc_index] = out_sc_pale[k].Std();
			}
	            if(out_sc_yell[k].HowMany() > 0 && CellTypes[k].sc_index >= 0) {
		        sc_yellmat[row][CellTypes[k].sc_index] = out_sc_yell[k].Mean()*double(out_sc_yell[k].HowMany())/NYELLOW;
                        sc_yell_sd[row][CellTypes[k].sc_index] = out_sc_yell[k].Std();
			}
		    }
	        }
	    }
	    
	in_row_ids.clear();   in_rows.clear();  in_rows_dat.clear();  in_row_sums.clear();  in_rows_recip.clear();
	out_row_ids.clear(); out_rows.clear(); out_rows_dat.clear(); out_row_sums.clear(); out_rows_recip.clear();
        in_row_area.clear(); out_row_area.clear(); in_row_vols.clear(); out_row_vols.clear();
        in_row_inc.clear(); out_row_inc.clear();
        out_stats.clear(); out_stats.resize(NumTypes, MeanStd());
         in_stats.clear();  in_stats.resize(NumTypes, MeanStd());
         out_7col.clear();  out_7col.resize(NumTypes, MeanStd());
         out_pale.clear();  out_pale.resize(NumTypes, MeanStd());
         out_yell.clear();  out_yell.resize(NumTypes, MeanStd());
         out_home.clear();  out_home.resize(NumTypes, MeanStd());
         out_hh  .clear();  out_hh  .resize(NumTypes, MeanStd());
          in_7col.clear();   in_7col.resize(NumTypes, MeanStd());
        in_pale.clear();   in_pale.resize(NumTypes, MeanStd());
        in_yell.clear();   in_yell.resize(NumTypes, MeanStd());
        in_sc_pale.clear();   in_sc_pale.resize(NumTypes, MeanStd());
        in_sc_yell.clear();   in_sc_yell.resize(NumTypes, MeanStd());
        out_sc_pale.clear();   out_sc_pale.resize(NumTypes, MeanStd());
        out_sc_yell.clear();   out_sc_yell.resize(NumTypes, MeanStd());
        for(int j=0; j<8; j++) 
	    for(int k=0; k<8; k++)
		in_cross_col[j][k] = out_cross_col[j][k] = 0;
	printf("\n");
	}
    if (ii == M)
	break;
    printf("----------------- Starting analysis of '%s'\n", names[i]);
    cc = color(names[i]);
    printf("----------------- Color is %d %x\n", cc,cc);

    int sos = 0;  // sum of strengths
    for(int j=0; j<M; j++) {
	sa[j].index = j;
        int str = from_to[i][j];
	sa[j].val = str;
	sos += str;
        if (str > 0) {
            int tii = TypeIndex(names[i]);
            int tij = TypeIndex(names[j]);
            if (tij == 13)
		printf("Connection to type 13 (Mi1) of strength %d\n", str);
	    out_stats[TypeIndex(names[j])].Element(str);
            if (IsPale(names[j]) && IsPale(names[i]) && IsColumnar(names[j]) && IsColumnar(names[i]))    // Pale to Pale
	        out_pale[TypeIndex(names[j])].Element(str);
            if (IsYell(names[j]) && IsYell(names[i]) && IsColumnar(names[j]) && IsColumnar(names[i]))    // Yellow to Yellow
	        out_yell[TypeIndex(names[j])].Element(str);
            if (IsHome(names[i]) )                                                                       // Home to any other
	        out_home[TypeIndex(names[j])].Element(str);
            if (IsHome(names[i])  && IsHome(names[j]) )                                                  // Home to home
	        out_hh[TypeIndex(names[j])].Element(str);
            if (Is7col(names[j]) && Is7col(names[i]) && IsColumnar(names[j]) && IsColumnar(names[i]))    // within 7 to within 7
	        out_7col[TypeIndex(names[j])].Element(str);
            if (Column(names[j]) != ' ' && Column(names[i]) == Column(names[j]) && IsPale(names[j])) {     // in a single pale column
                printf("Adding pale %s -> %s, strength %d, col '%c'\n", names[i], names[j], str, Column(names[j]) );
                out_sc_pale[TypeIndex(names[j])].Element(str);
		}
            if (Column(names[j]) != ' ' && Column(names[i]) == Column(names[j]) && IsYell(names[j])) {    // in a single yellow column
                printf("Adding yell %s -> %s, strength %d, col '%c'\n", names[i], names[j], str, Column(names[j]) );
                out_sc_yell[TypeIndex(names[j])].Element(str);
                }
            out_cross_col[NumColumn(names[i])][NumColumn(names[j])] += str;
            }
	}
    sort(sa.begin(), sa.end());
    if (sa[0].val > 0) {
        printf("cell %11s -> ", names[i]);
        out_row_ids.push_back(i);
        out_row_sums.push_back(sos);
        out_row_vols.push_back(volumes[i]);
        out_row_area.push_back(  areas[i]);
        out_row_inc .push_back(incomplete[i]);
        vector<int> one_row;
        vector<int> one_row_dat;
        vector<int> one_row_recip;
        for(int k=0; k<NUM_SHOWN && sa[k].val > 0; k++) {
	    int ind = sa[k].index;
            one_row.push_back(ind);
            one_row_dat.  push_back(from_to[i][ind]);
            one_row_recip.push_back(from_to[ind][i]);
	    }
        out_rows.push_back(one_row);
        out_rows_dat.push_back(one_row_dat);
        out_rows_recip.push_back(one_row_recip);
        printf("\n");
	}
    // now do the same for incoming
    sos = 0;
    for(int j=0; j<M; j++) {
	sa[j].index = j;
        int str = from_to[j][i];
	sa[j].val = str;
	sos += str;
        if (str > 0) {
            int tii = TypeIndex(names[i]);
            int tij = TypeIndex(names[j]);
	    in_stats[TypeIndex(names[j])].Element(str);
            if (IsPale(names[j]) && IsPale(names[i]) && IsColumnar(names[j]) && IsColumnar(names[i]))
	        in_pale[TypeIndex(names[j])].Element(str);
            if (IsYell(names[j]) && IsYell(names[i]) && IsColumnar(names[j]) && IsColumnar(names[i]))
	        in_yell[TypeIndex(names[j])].Element(str);
            if (Is7col(names[j]) && Is7col(names[i]) && IsColumnar(names[j]) && IsColumnar(names[i]))
	        in_7col[TypeIndex(names[j])].Element(str);
            if (Column(names[j]) != ' ' && Column(names[i]) == Column(names[j]) && IsPale(names[j]))
                in_sc_pale[TypeIndex(names[j])].Element(str);
            if (Column(names[j]) != ' ' && Column(names[i]) == Column(names[j]) && IsYell(names[j]))
                in_sc_yell[TypeIndex(names[j])].Element(str);
            in_cross_col[NumColumn(names[i])][NumColumn(names[j])] += str;
            }
	}
    sort(sa.begin(), sa.end());
    if (sa[0].val > 0) {
        printf("cell %11s <- ", names[i]);
        in_row_ids.push_back(i);
        in_row_sums.push_back(sos);
        in_row_vols.push_back(volumes[i]);
        in_row_area.push_back(  areas[i]);
        in_row_inc .push_back(incomplete[i]);
        vector<int> one_row;
        vector<int> one_row_dat;
        vector<int> one_row_recip;
        for(int k=0; k<NUM_SHOWN && sa[k].val > 0; k++) {
	    int ind = sa[k].index;
            one_row.push_back(ind);
            one_row_dat.  push_back(from_to[ind][i]);
            one_row_recip.push_back(from_to[i][ind]);
	    }
        in_rows.push_back(one_row);
        in_rows_dat.push_back(one_row_dat);
        in_rows_recip.push_back(one_row_recip);
        printf("\n");
	}
    }
return 0;  // Stop at this point for MB
// Write the table equivalent to the single cell case.
int old[] = {
  0,119,  0,116,100, 52,  0,   0,  0, 61,  0,  0,  0,  0,  0,  0,  0,   0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
  8,  6,199, 36,  8,  6, 17,   0,  0,  2,  0,  0,  0,  7,  0,  0, 10,   0,  0,  5,  2,  0,  0,  1,  0,  5,  1,
  0,  0, 12,  0,  0,  0,  0,   0,  0,  4,  0,  0,  0,  1,  0,  0,  0,   0,  0,  0, 16,  0,  0,  0,  0,  0,  0,
  7,  1, 86,  7,  0,  2,  0,   0,  0,  0,  0,  0,  0,  5,  0,  6,  1,   0,  0,  0,  0,  0,  0,  0,  0,  0,  1,
 20, 37,  0, 39,  3, 12,  0,   4,  0, 12,  0,  0,  0, 14,  0,  0,  5,   0,  0, 33,  0,  0,  0,  0,  0, 11,  2,
 44, 22,  0,  0, 36,  0,  0,  26,  8,  0,  0,  1,  0,  0,  8,  0,  0,   0, 15,  0,  0,  0,  0,  0,  0,  0,  7,
  0,  2,  0,  2,  0,  0,  0,   0,  0,  2,  0,  0,  0,  0,  0,  7,  0,   1,  0,  0,  0,  0,  1,  0,  0,  0,  0,

  0,  0,  0,  0, 35,  0,  0,   0,153,  0,106, 56, 67,  0,  0,  0,  0,   0,  0,  0,  0,  0, 10,  0,  0,  0,  0,
  0,  1,  0,  0, 12,  0,  5,   0,  0,  0,  1,  0,  3,  1,  0,  0,  6,   0,  0,  0,  4,  0,  2,  0,  1,  0,  0,
  0,  0,  5,  0,  9,  0,  0,  37, 18,  0,  0, 62, 20,  3,  0,  0,  0,   0,  0,  0,  7,  0, 10,  0,  3,  0,  0,
  0,  0,  0,  0, 12,  0,  1,   0,  0,  4,  0,  0,  3, 10,  1,  1,  2,   0,  0,  1,  6,  0,  0,  0,  1,  0,  0,
  0,  0,  0,  0,  0,  0,  0,   0,  0,  0,  0,  0,  0,  0,  0,  0,  0,   0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
  0,  1,  0,  1,  0,  2,  2,   0,  1,  2,  1,  0,  0,  0,  0, 13,  0,   0,  0,  2,  0,  0,  0,  0,  0,  0,  0,
  0,  5,  0,  0,  0,  0,  0,   0,  0,  0,  0,  0,  1,  1,  0,  4,  2,   0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
  0,  0,  0,  0,  0,  0,  0,   0,  0,  0, 15,  0,  5,  0,  1,  0,  0,   0,  0,  0,  2,  0,  0,  0,  0,  0,  0,
  0,  0,  0,  0,  0,  0,  0,   0,  0,  0,  0,  0,  2,  0,  0,  0,  1,   0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
  0,  0,  0,  2,  0,  0,  0,   0,  1,  0,  0,  0,  0,  0,  0,  0,  0,   0,  0,  0,  0,  0,  0,  0,  0,  0,  0,

  6, 20,  0,  0,  0,  0,  0,   0,  0,  0,  0,  0,  0,  0,  0,  0,  0,   0,  0, 25,  1, 11, 11,  0,  0,  9, 11,
  0, 14,  0,  0,  0,  0,  0,   0,  0,  9,  0,  0,  0,  0,  0,  0,  0,   0,  0,  0, 30,  0, 20,  0, 17,  0,  0,
  0,  0,  4,  1,  1,  0,  0,   0,  0,  0,  0,  0,  0,  1,  0,  2,  0,   0,  0,  0, 18,  0,  7,  0,  7,  0,  0,
  0,  0,  6,  0,  0,  0,  0,   0,  5,  0,  1,  0,  2,  3,  0,  1,  4,   0,  0,  5,  0,  0,  0,  0,  0,  0,  0,
  0,  0,  0,  0,  0,  0,  0,   3,  0,  0,  0,  0,  0,  0,  0,  0,  0,   5,  3,  0,  1,  0,  0, 35,  0,  1,  0,
  0,  0,  0,  0,  0,  0,  1,   0,  0,  0,  0,  1,  0,  0,  0,  0,  0,   0,  0,  0,  5,  0,  1,  0,  3,  0,  0,
  0,  0,  0,  0,  0,  0,  0,   0,  0,  0,  0,  0,  1,  0,  0,  0,  0,   0,  0,  1,  0,  0,  0,  0,  0,  3,  0,
  0,  0,  0,  0,  0,  0,  0,   0,  0,  0,  0,  0,  0,  0,  0,  0,  0,   0,  0,  0,  0,  0,  0,  0,  0,  1,  0,
  0,  0,  0,  0,  0,  0,  0,   0,  0,  0,  0,  0,  0,  0,  0,  0,  0,   0,  0, 11,  0,  0,  0,  0,  0,  2,  0,
  0,  0,  0,  0,  0,  0,  0,   0,  0,  0,  0,  0,  3,  0,  0,  1,  0,   0,  1,  6,  0,  0,  0,  0,  0,  1,  2,
};

// make an array that tells the cell type for each row column;
//
int old_index[NS];
for(int i=0; i<NumTypes; i++) {
    for(int k=0; k<NumTypes; k++)
	if (CellTypes[k].sc_index == i)
	    old_index[i] = k;
    }

for(int row=0; row<NS; row++) {
    for(int col=0; col<NS; col++)
	printf("%4d", int(scmat[row][col]+0.5) );
    printf("\n");
    }
FILE *fp = fopen("table.html","w");
FILE *f2 = fopen("table.csv", "w");
if (fp == NULL || f2 == NULL) {
    printf("Could not open output files.\n");
    return 42;
    }
    
fprintf(fp, "<!DOCTYPE html PUBLIC \"-//W3C//DTD HTML 4.01 Transitional//EN\">\n");
fprintf(fp, "<html>\n");
fprintf(fp, "  <head>\n");
fprintf(fp, "    <meta http-equiv=\"content-type\" content=\"text/html;\n");
fprintf(fp, "      charset=windows-1252\">\n");
fprintf(fp, "    <title> Table </title>\n");
fprintf(fp, "  </head>\n");
fprintf(fp, "  <body>\n");
// Contents go here
fprintf(fp, "<p>Comparison of new results vs Takemura et al., 2013. (bottom).</p>"
            "<p>The top number is the number of connections from the home column to all columns.</p>"
            "<p>The second number is from the home column to the home column only.</p>"
            "<p>The bottom number is the strength reported in Takemura, et. al., in 2013. </p>\n");
int blue = 0xFF;
fprintf(fp, "<TABLE BORDER=4 CELLSPACING=4 CELLPADDING=4>\n");
fprintf(fp, "<TR><TD> ");
for(int col=0; col<NS; col++) {
    fprintf(fp, "<TD>%s\n", CellTypes[old_index[col]].name);
    fprintf(f2,",%s", CellTypes[old_index[col]].name);
    }
fprintf(f2,"\n");    
fprintf(fp, "</TD></TR>\n");

int sum_new = 0;
int sum_old = 0;
for(int row=0; row<NS; row++) {
    int red = 0xFF0000;
    if (row <= 16) red = 0xBB0000;
    if (row <= 6) red = 0xDD0000;
    fprintf(fp, "<TR><TD>%s", CellTypes[old_index[row]].name);
    fprintf(f2, "%s", CellTypes[old_index[row]].name);
    for(int col=0; col<NS; col++) {
        int green = 0xFF00;
        if (col <= 16) green = 0xBB00;
        if (col <= 6) green = 0xDD00;
        double mult = (row == 2 || col == 2) ? 1.0 : 1.0;  // Multiply row 2 (T4) by 4, since ours is per cell and old paper lumped.
	int newv = int(homemat[row][col]*mult + 0.5);
        int hhv  = int(hhmat[row][col] + 0.5);             // value only to-from home column
        int oldv = old[row*NS+col];
        bool interesting = (oldv == 0 && newv >= 3) || (oldv >= 3 && newv == 0);
        int avg = (newv + oldv) /2;
        interesting = interesting || ((avg >= 5) && (newv > oldv*3 || oldv > 3*newv));
        const char *font = (interesting && newv < oldv) ? "<FONT COLOR=\"#FF0000\">" : ""; // red if new value is less
	sum_new += newv;
        sum_old += oldv;
	fprintf(fp, "<TD BGCOLOR=\"#%06x\">%s%4d<br>%4d<br>%4d", red+green+(interesting?0:blue), font, newv, hhv, oldv);
        fprintf(f2,",%d", newv);
	}
    fprintf(fp, "</TD></TR>\n");
    fprintf(f2, "\n");
    }
fprintf(fp, "</TABLE>\n");
fprintf(fp, "<p>Total connections in old table %d, new table %d, ratio %.3f</p>\n", sum_old, sum_new, double(sum_new)/sum_old );
fprintf(fp, "</body>\n");
fprintf(fp, "</html>\n");
fclose(fp);
fclose(f2);

//---------------------------- Now write differences between home, pale, and yellow columns
//
fp = fopen("table2.html","w");
f2 = fopen("table2.csv", "w");
if (fp == NULL || f2 == NULL) {
    printf("Could not open output files.\n");
    return 42;
    }
    
fprintf(fp, "<!DOCTYPE html PUBLIC \"-//W3C//DTD HTML 4.01 Transitional//EN\">\n");
fprintf(fp, "<html>\n");
fprintf(fp, "  <head>\n");
fprintf(fp, "    <meta http-equiv=\"content-type\" content=\"text/html;\n");
fprintf(fp, "      charset=windows-1252\">\n");
fprintf(fp, "    <title> Table </title>\n");
fprintf(fp, "  </head>\n");
fprintf(fp, "  <body>\n");
// Contents go here
fprintf(fp,"<p> This table shows four numbers per box.  The top left is over all identified cells (including fragments).  The top right "
"is all the connections within our seven column set.  The bottem left is within the pale columns, and the bottom right within the"
" yellow columns.  Yellow boxes are when the average count is >= 8 and the difference at least 1.3</p>\n");
fprintf(fp, "<TABLE BORDER=4 CELLSPACING=4 CELLPADDING=4>\n");
fprintf(fp, "<TR><TD> ");
for(int col=0; col<NS; col++) {
    fprintf(fp, "<TD>%s\n", CellTypes[old_index[col]].name);
    fprintf(f2,",%s", CellTypes[old_index[col]].name);
    }
fprintf(f2,"\n");    
fprintf(fp, "</TD></TR>\n");

for(int row=0; row<NS; row++) {
    int red = 0xFF0000;
    if (row <= 16) red = 0xBB0000;
    if (row <= 6) red = 0xDD0000;
    fprintf(fp, "<TR><TD>%s", CellTypes[old_index[row]].name);
    fprintf(f2, "%s", CellTypes[old_index[row]].name);
    for(int col=0; col<NS; col++) {
        int green = 0xFF00;
        if (col <= 16) green = 0xBB00;
        if (col <= 6) green = 0xDD00;
        double mult = (row == 2 || col == 2) ? 4.0 : 1.0;  // Multiply row 2 (T4) by 4, since ours is per cell and old paper lumped.
	int newv = int(scmat[row][col]*mult + 0.5);
	int new7 = int(cen7mat[row][col]*mult + 0.5);
	int newp = int(palemat[row][col]*mult + 0.5);
	int newy = int(yellmat[row][col]*mult + 0.5);
        int avg = (newp + newy)/2;
        bool interesting = avg >= 8 && (newp > 1.3*newy || newy > 1.3*newp);
        int oldv = 0;
        const char *font = (interesting && newy > newp) ? "<FONT COLOR=\"#FF0000\">" : ""; // red if yellow is greater
        if (max(newv, max(new7, max(newp,newy))) < 100)
	    fprintf(fp, 
            "<TD BGCOLOR=\"#%06x\">%s<pre>%2d%3d<br>%2d%3d</pre>", red+green+(interesting?0:blue), font, newv, new7, newp, newy);
        else
	    fprintf(fp,
            "<TD BGCOLOR=\"#%06x\">%s<pre>%3d%4d<br>%3d%4d</pre>", red+green+(interesting?0:blue), font, newv, new7, newp, newy);
        fprintf(f2,",%d", newv);
	}
    fprintf(fp, "</TD></TR>\n");
    fprintf(f2, "\n");
    }
fprintf(fp, "</TABLE>\n");
fprintf(fp, "</body>\n");
fprintf(fp, "</html>\n");
fclose(fp);
fclose(f2);

// --------------------------------   Write out the data to/from exact same column
fp = fopen("table3.html","w");
f2 = fopen("table3.csv", "w");
if (fp == NULL || f2 == NULL) {
    printf("Could not open output files.\n");
    return 42;
    }
    
fprintf(fp, "<!DOCTYPE html PUBLIC \"-//W3C//DTD HTML 4.01 Transitional//EN\">\n");
fprintf(fp, "<html>\n");
fprintf(fp, "  <head>\n");
fprintf(fp, "    <meta http-equiv=\"content-type\" content=\"text/html;\n");
fprintf(fp, "      charset=windows-1252\">\n");
fprintf(fp, "    <title> Table </title>\n");
fprintf(fp, "  </head>\n");
fprintf(fp, "  <body>\n");
// Contents go here
fprintf(fp, "<p>This table has 6 entries per box.  The top two are the total entries per type (including partial cells) and"
" then connections only between cells in the main 7 columns.  The second row has the average numbers within the pale (first)"
" and yellow (second) numbers.  The third row has the standard deviation of the pale and yellow columns.  This is to help"
" determine if the difference is significant.  Yellow squares are where the difference is more than two (averaged) standard"
" deviations.</p>\n");
fprintf(fp, "<TABLE BORDER=4 CELLSPACING=4 CELLPADDING=4>\n");
fprintf(fp, "<TR><TD> ");
for(int col=0; col<NS; col++) {
    fprintf(fp, "<TD>%s\n", CellTypes[old_index[col]].name);
    fprintf(f2,",%s", CellTypes[old_index[col]].name);
    }
fprintf(f2,"\n");    
fprintf(fp, "</TD></TR>\n");

for(int row=0; row<NS; row++) {
    int red = 0xFF0000;
    if (row <= 16) red = 0xBB0000;
    if (row <= 6) red = 0xDD0000;
    fprintf(fp, "<TR><TD>%s", CellTypes[old_index[row]].name);
    fprintf(f2, "%s", CellTypes[old_index[row]].name);
    for(int col=0; col<NS; col++) {
        int green = 0xFF00;
        if (col <= 16) green = 0xBB00;
        if (col <= 6) green = 0xDD00;
        double mult = (row == 2 || col == 2) ? 4.0 : 1.0;  // Multiply row 2 (T4) by 4, since ours is per cell and old paper lumped.
	int newv = int(scmat[row][col]*mult + 0.5);
	int new7 = int(cen7mat[row][col]*mult + 0.5);
	int newp = int(sc_palemat[row][col]*mult + 0.5);
	int newy = int(sc_yellmat[row][col]*mult + 0.5);
        int avg = (newp + newy)/2;
        double std_avg = (sc_pale_sd[row][col] + sc_yell_sd[row][col])/2.0;
        bool interesting = avg >= 8 && abs(newp-newy) >= 2*std_avg;
        const char *font = (interesting && newy > newp) ? "<FONT COLOR=\"#FF0000\">" : ""; // red if yellow is greater
        if (max(newv, max(new7, max(newp,newy))) < 100)
	    fprintf(fp, 
            "<TD BGCOLOR=\"#%06x\">%s<pre>%2d%3d<br>%2d%3d<br>%.1f %.1f</pre>", red+green+(interesting?0:blue), font, 
            newv, new7, newp, newy, sc_pale_sd[row][col], sc_yell_sd[row][col]);
        else
	    fprintf(fp,
            "<TD BGCOLOR=\"#%06x\">%s<pre>%3d%4d<br>%3d%4d<br>%.1f %.1f</pre>", red+green+(interesting?0:blue), font, 
            newv, new7, newp, newy, sc_pale_sd[row][col], sc_yell_sd[row][col]);
        fprintf(f2,",%d", newv);
	}
    fprintf(fp, "</TD></TR>\n");
    fprintf(f2, "\n");
    }

fprintf(fp, "<TR><TD> ");
for(int col=0; col<NS; col++) {
    fprintf(fp, "<TD>%s\n", CellTypes[old_index[col]].name);
    }
fprintf(fp, "</TD></TR>\n");

fprintf(fp, "</TABLE>\n");
fprintf(fp, "</body>\n");
fprintf(fp, "</html>\n");
fclose(fp);
fclose(f2);
// ----------------------------- Write tables of cell types sorted by deviation
sa.resize(CellStats.size());  // make the sort array the right size
fp = fopen("table4.html", "w");
if (fp == NULL) {
    printf("Cannot open table4.html\n");
    return 42;
    }
fprintf(fp, "<!DOCTYPE html PUBLIC \"-//W3C//DTD HTML 4.01 Transitional//EN\">\n");
fprintf(fp, "<html>\n");
fprintf(fp, "  <head>\n");
fprintf(fp, "    <meta http-equiv=\"content-type\" content=\"text/html;\n");
fprintf(fp, "      charset=windows-1252\">\n");
fprintf(fp, "    <title> Table </title>\n");
fprintf(fp, "  </head>\n");
fprintf(fp, "  <body>\n");
fprintf(fp, "<p>These tables are sorted by stereotypy, with the most stereotyped first.  The number reported is the standard deviation divided by "
" the mean.  The cells are color coded according to the volume stereotype, divided into three bands.  'Total' is the total numbers of"
" each cell type we identified.  However, only the cells that are complete (not cut off at the edge) are used in the calculation of"
" the deviation.</p>\n");
// first by volume
printf("Sorted by volume stereotypy\n");
for(int i=0; i<sa.size(); i++) {
    sa[i].index = i;
    sa[i].val = CellStats[i].volume_dev;
    };
sort(sa.begin(), sa.end());
// set the colors.  First count the ones eligible to be printed.
int eligible = 0;
for(int i=0; i<CellStats.size(); i++) {
    if (CellStats[i].used_stats >= 2)
	eligible++;
    }
int cnt = 0;
for(int i=0; i<sa.size(); i++) {
    int j = sa[i].index;
    if (CellStats[j].how_many < 2)
	continue; // Not eligible
    cnt++;
    if (cnt < eligible/3)
	CellStats[j].color = 0xCC0000; // red
    else if (cnt < (eligible*2)/3 )
	CellStats[j].color = 0xCCCC00; // yellow
    else
	CellStats[j].color = 0x00CC00; // green
    }
PrintStereo(fp, CellStats, sa, "Sorted by variance of volume");
printf("Sorted by input stereotypy\n");
for(int i=0; i<sa.size(); i++) {
    sa[i].index = i;
    sa[i].val = CellStats[i].in_dev;
    };
sort(sa.begin(), sa.end());
PrintStereo(fp, CellStats, sa, "Sorted by variance of number of inputs");
printf("Sorted by output stereotypy\n");
for(int i=0; i<sa.size(); i++) {
    sa[i].index = i;
    sa[i].val = CellStats[i].out_dev;
    }
sort(sa.begin(), sa.end());
PrintStereo(fp, CellStats, sa, "Sorted by variance of number of outputs");
printf("Sorted by surface area stereotypy\n");
for(int i=0; i<sa.size(); i++) {
    sa[i].index = i;
    sa[i].val = CellStats[i].area_dev;
    }
sort(sa.begin(), sa.end());
PrintStereo(fp, CellStats, sa, "Sorted by variance of surface area");
fprintf(fp, "</body>\n");
fprintf(fp, "</html>\n");
fclose(fp);

//------------------------------- Analyze the correlation between rows and columns
//We will write two files, one for inputs and one for outputs.  Each will have two table - one with all cells (including
//fragments) and one with cells known to be in the 7 columns.
fp = fopen("InputCorr.html", "w");
if (fp == NULL) {
    printf("Cannot open InputCorr.html\n");
    return 42;
    }
fprintf(fp, "<!DOCTYPE html PUBLIC \"-//W3C//DTD HTML 4.01 Transitional//EN\">\n");
fprintf(fp, "<html>\n");
fprintf(fp, "  <head>\n");
fprintf(fp, "    <meta http-equiv=\"content-type\" content=\"text/html;\n");
fprintf(fp, "      charset=windows-1252\">\n");
fprintf(fp, "    <title> Table </title>\n");
fprintf(fp, "  </head>\n");
fprintf(fp, "  <body>\n");
fprintf(fp, "<p>This table gives the normalized inner product between the vector of inputs to each of the modular neuron type."
"   If they both got all their inputs, in the same proportion, from exactly the same neurons, the inner product is 1.  If the "
"inputs come from completely different sources, the inner product is zero.  This table includes all cells, including those in "
"outside the home 7, and partially reconstructed cells.  This makes the data noisier but allows non-columnar cells such as "
"Tm3 and Tm4 to be included.</p>\n");


FILE *fo = fopen("OutputCorr.html", "w");
if (fp == NULL) {
    printf("Cannot open OutputCorr.html\n");
    return 42;
    }
fprintf(fo, "<!DOCTYPE html PUBLIC \"-//W3C//DTD HTML 4.01 Transitional//EN\">\n");
fprintf(fo, "<html>\n");
fprintf(fo, "  <head>\n");
fprintf(fo, "    <meta http-equiv=\"content-type\" content=\"text/html;\n");
fprintf(fo, "      charset=windows-1252\">\n");
fprintf(fo, "    <title> Table </title>\n");
fprintf(fo, "  </head>\n");
fprintf(fo, "  <body>\n");
fprintf(fo, "<p>This table gives the normalized inner product between the vector of outputs from each of the modular neuron type."
"   If they both send all outputs, in the same proportion, to exactly the same neurons, the inner product is 1.  If the "
"go to completely different neurons, the inner product is zero.  This table includes all cells, including those in "
"outside the home 7, and partially reconstructed cells.  This makes the data noisier but allows non-columnar cells such as "
"Tm3 and Tm4 to be included.</p>\n");
#define MATRIX scmat
fprintf(fo, "<TABLE BORDER=4 CELLSPACING=4 CELLPADDING=4>\n");
fprintf(fo, "<TR><TD> ");
for(int col=0; col<NS; col++) {
    fprintf(fo, "<TD>%s\n", CellTypes[old_index[col]].name);
    }
fprintf(fo, "</TD></TR>\n");
for(int r1 = 0; r1 < MATRIX.size(); r1++) {
    fprintf(fo, "<TR><TD>%s", CellTypes[old_index[r1]].name);
    for(int r2 = 0; r2 < MATRIX.size(); r2++) {
	//printf(" coeff= %d %d = %5.2f\n",r1, r2, Inner(MATRIX[r1], MATRIX[r2]) );
	double ip = Inner(MATRIX[r1], MATRIX[r2]);
	int color = (ip > 0.695 && r1 != r2) ? 0xAAFFAA : ((ip > 0.595 && r1 != r2)? 0xFFFFAA : 0xFFFFFF);
	if (r2 < r1)
	    fprintf(fo, "<TD>     ");
	else
	    fprintf(fo, "<TD BGCOLOR=\"#%06x\">%5.2f", color, Inner(MATRIX[r1], MATRIX[r2]) );
        }
    fprintf(fo,"</TD></TR>\n");
    }
fprintf(fo, "</TABLE>\n");

// now transpose the matrix
for(int r1 = 0; r1 < MATRIX.size(); r1++) {
    for(int r2 = 0; r2 < MATRIX.size(); r2++) {
        if (r1 < r2) {
	    double swap = MATRIX[r1][r2];
	    MATRIX[r1][r2] = MATRIX[r2][r1];
	    MATRIX[r2][r1] = swap;
	    }
	}
    }

// And print again
fprintf(fp, "<TABLE BORDER=4 CELLSPACING=4 CELLPADDING=4>\n");
fprintf(fp, "<TR><TD> ");
for(int col=0; col<NS; col++) 
    fprintf(fp, "<TD>%s\n", CellTypes[old_index[col]].name);
fprintf(fp, "</TD></TR>\n");
for(int r1 = 0; r1 < MATRIX.size(); r1++) {
    fprintf(fp, "<TR><TD>%s", CellTypes[old_index[r1]].name);
    for(int r2 = 0; r2 < MATRIX.size(); r2++) {
	//printf(" coeff= %d %d = %5.2f\n",r1, r2, Inner(MATRIX[r1], MATRIX[r2]) );
	double ip = Inner(MATRIX[r1], MATRIX[r2]);
	int color = (ip > 0.895 && r1 != r2) ? 0xAAFFAA : ((ip > 0.795 && r1 != r2)? 0xFFFFAA : 0xFFFFFF);
	if (r2 < r1)
	    fprintf(fp, "<TD>     ");
	else
	    fprintf(fp, "<TD BGCOLOR=\"#%06x\">%5.2f", color, Inner(MATRIX[r1], MATRIX[r2]) );
        }
    fprintf(fp,"</TD></TR>\n");
    }
fprintf(fp, "</TABLE>\n");
#undef MATRIX
#define MATRIX cen7mat
fprintf(fo, "<TABLE BORDER=4 CELLSPACING=4 CELLPADDING=4>\n");
fprintf(fo, "<TR><TD> ");
for(int col=0; col<NS; col++) {
    fprintf(fo, "<TD>%s\n", CellTypes[old_index[col]].name);
    }
fprintf(fo, "</TD></TR>\n");
for(int r1 = 0; r1 < MATRIX.size(); r1++) {
    fprintf(fo, "<TR><TD>%s", CellTypes[old_index[r1]].name);
    for(int r2 = 0; r2 < MATRIX.size(); r2++) {
	//printf(" coeff= %d %d = %5.2f\n",r1, r2, Inner(MATRIX[r1], MATRIX[r2]) );
	double ip = Inner(MATRIX[r1], MATRIX[r2]);
	int color = (ip > 0.695 && r1 != r2) ? 0xAAFFAA : ((ip > 0.595 && r1 != r2)? 0xFFFFAA : 0xFFFFFF);
	if (r2 < r1)
	    fprintf(fo, "<TD>     ");
	else
	    fprintf(fo, "<TD BGCOLOR=\"#%06x\">%5.2f", color, Inner(MATRIX[r1], MATRIX[r2]) );
        }
    fprintf(fo,"</TD></TR>\n");
    }
fprintf(fo, "</TABLE>\n");

// now transpose the matrix
for(int r1 = 0; r1 < MATRIX.size(); r1++) {
    for(int r2 = 0; r2 < MATRIX.size(); r2++) {
        if (r1 < r2) {
	    double swap = MATRIX[r1][r2];
	    MATRIX[r1][r2] = MATRIX[r2][r1];
	    MATRIX[r2][r1] = swap;
	    }
	}
    }

// And print again
fprintf(fp, "<TABLE BORDER=4 CELLSPACING=4 CELLPADDING=4>\n");
fprintf(fp, "<TR><TD> ");
for(int col=0; col<NS; col++) {
    fprintf(fp, "<TD>%s\n", CellTypes[old_index[col]].name);
    }
fprintf(fp, "</TD></TR>\n");
for(int r1 = 0; r1 < MATRIX.size(); r1++) {
    fprintf(fp, "<TR><TD>%s", CellTypes[old_index[r1]].name);
    for(int r2 = 0; r2 < MATRIX.size(); r2++) {
	//printf(" coeff= %d %d = %5.2f\n",r1, r2, Inner(MATRIX[r1], MATRIX[r2]) );
	double ip = Inner(MATRIX[r1], MATRIX[r2]);
	int color = (ip > 0.895 && r1 != r2) ? 0xAAFFAA : ((ip > 0.795 && r1 != r2)? 0xFFFFAA : 0xFFFFFF);
	if (r2 < r1)
	    fprintf(fp, "<TD>     ");
	else
	    fprintf(fp, "<TD BGCOLOR=\"#%06x\">%5.2f", color, Inner(MATRIX[r1], MATRIX[r2]) );
        }
    fprintf(fp,"</TD></TR>\n");
    }
fprintf(fp, "</TABLE>\n");

fprintf(fp, "</body>\n"); fprintf(fo, "</body>\n");
fprintf(fp, "</body>\n"); fprintf(fo, "</html>\n");
fclose(fp);               fclose(fo);

// Find all the names that occur in all 7 columns
vector<int> counts(NumTypes, 0);
for(int i=0; i<names.size(); i++) {
    if (isdigit(names[i][0]))
	continue;
    if (Is7col(names[i])) {
        int ty = TypeIndex(names[i]);
        counts[ty]++;
        printf("Cell %s is part of 7 main columns; count now %d\n", names[i], counts[ty]);
        }
    }
int Nrac = 0;  // number represented in all columns, and their names
vector<const char *> all_col_names;
for(int i=0; i<NumTypes; i++) {
    // Here we have two options, either 6/7 or all 7 required
    //if ((counts[i] == 7 || counts[i] == 6) && CellTypes[i].new_col) {
    if (counts[i] == 7 && CellTypes[i].new_col) {
	printf("Type %s occurs in all 7 columns\n", CellTypes[i].name);
        counts[i] = Nrac;
        all_col_names.push_back(CellTypes[i].name);
        Nrac++;
	}
    else
        counts[i] = -1;
    }
printf("There are %d types that occur in every column\n", Nrac);
// At this point, counts[i] is the index of cell type i, in a table of celltypes.
// Construct 'ic' where ic will be a Nrac*7 matrix, looking like
//     L1H L1a L1b L1c L1d L1e L1f C2H c2a C2b C2c C2d C2e C2f ....
// L1H  s1  s2  s3
// L1a
// L1b
// L1c
// ...
// where the 's' entries are strengths.
int grand_total = 0;
vector<int> subset_index(Nrac*7, -1);
vector< vector<int> > ic(Nrac*7, vector<int> (Nrac*7, 0));
for(int s=0; s<names.size(); s++) {
    int ty = TypeIndex(names[s]);
    if (counts[ty] < 0 || !Is7col(names[s]) )
	continue;
    int is = counts[ty];
    int cs = NumColumn(names[s]) - 1; // -1 to convert to 0 based.
    subset_index[is*7+cs] = s;
    for(int t=0; t<names.size(); t++) {
	if (from_to[s][t] > 0) {
	    int ty2 = TypeIndex(names[t]);
	    if (counts[ty2] < 0 || !Is7col(names[t]) )
		continue;
            int it = counts[ty2];
            int ct = NumColumn(names[t]) - 1;
	    printf("From %s to %s, weight %d %6d %3d %3d %3d\n", names[s], names[t], from_to[s][t], is, cs, it, ct );
            if (is < 0 || is >= Nrac || cs < 0 || cs > 6 || it < 0 || it >= Nrac || ct < 0 || ct > 6) {
		printf("Bogons!  %d %d %d %d\n", is, cs, it, ct);
		return 42;
		}
            ic[is*7+cs][it*7+ct] = from_to[s][t];
            grand_total += from_to[s][t];
	    }
	}
    }
// Write the core connectome as a comma separated file
printf(" ,");
for(int i=0; i<Nrac*7; i++)
    printf("%s%s", names[subset_index[i]], i != Nrac*7-1 ? "," : "");
printf("\n");
for(int from = 0; from<7*Nrac; from++) {
    printf("%s,", names[subset_index[from]]);
    for(int to=0; to<Nrac*7; to++)
	printf("%d%s", ic[from][to], to != Nrac*7-1 ? "," : "");
    printf("\n");
    }
// Variability of TBar counts
for(int ml=0; ml<=10; ml++) {
    char name[128];
    sprintf(name,"/tmp/tb.%d", ml);
    FILE *f = fopen(name, "w");
    if (f == NULL) {
	printf("Could not open %s\n", name);
	return 42;
	}
    for(int i=0; i<Nrac; i++) {
	MeanStd m;
	for(int j=0; j<7; j++) {
	    int id = subset_index[7*i+j];
	    printf("Tbar: %s %d\n", names[id], Tbar_count_by_layer[id][ml]);
	    m.Element(Tbar_count_by_layer[id][ml]);
	    }
	printf("Tbar:: %.2f %.2f\n", m.Mean(), m.Std()*m.Std() );
	if (m.Sum() >= 1)
            fprintf(f, "%f %f \"%s\"\n", m.Mean(), m.Std()*m.Std(), BaseName(names[subset_index[7*i]]) );
	}
    fclose(f);
    }
printf("Grand total of connections among complete subset is %d\n", grand_total);
// read CoreStats.txt if available.
vector<double>      CoreVolume(Nrac*7, 0.0);
vector<double> CoreSurfaceArea(Nrac*7, 0.0);
vector<vector<double> >CoreOverlaps(Nrac*7, vector<double>(Nrac*7,0.0));
FILE *fover = fopen("CoreStats.txt", "r");
if (fover != NULL) {
    int N;
    fscanf(fover, "%d", &N);
    if (N != Nrac*7) {
	printf("Size mismatch %d %d\n", N, Nrac*7);
	return 42;
	}
    for(int y=0; y<N; y++) {
        int bid;
	if (fscanf(fover, "%d %lf %lf", &bid, &CoreVolume[y], &CoreSurfaceArea[y] ) != 3) {
	    printf("CoreStats error #1\n");
	    return 42;
	    }
        if (bid != rev_map[subset_index[y]]) {
	    printf("Mismatch in CoreStats.txt . Expected body %d, got %d\n", rev_map[subset_index[y]], bid);
	    return 42;
	    }
	}
    for(int y=0; y<N; y++) {
	for(int x=0; x<N; x++) {
	    if (fscanf(fover, "%lf", &(CoreOverlaps[y][x]) ) != 1) {
	        printf("CoreStats error #2\n");
		return 42;
	        }
	    }
	}
    fclose(fover);
    }
// Compute completeness for core connectome cells only
printf("Completeness by PSDs identified for core connectome\n");
MeanStd core_comp;
MeanStd core_comp_w;
double sum_of_weights = 0.0;
for(int i=0; i<subset_index.size(); i++) {
    int from = subset_index[i];
    if (isdigit(names[from][0]))
	return 42;  // should never happen
    int known = 0; int unknown = 0;
    for(int to=0; to<names.size(); to++) {
	int str = from_to[from][to];
	if (str == 0) 
	    continue;
	if (isdigit(names[to][0]))
	    unknown += str;
	else
	    known += str;
	}
    double pct = known/double(known+unknown);
    printf(" Cell %15s: %4d known, %4d unknown, %.2f%% complete\n", names[from], known, unknown, pct*100.0);
    core_comp.Element(pct);
    core_comp_w.Element(pct*(known+unknown));
    sum_of_weights += (known + unknown);
    }
printf("Core completeness %f raw, %f weighted\n", core_comp.Mean(), core_comp_w.Sum()/sum_of_weights);
printf("End of Completeness\n");
// Now read the same thing for all 10 medulla layers.  entry [0] will not be used.
const int Mlayers = 10;
vector<vector<double> >      VolumeByLayer(Mlayers+1, vector<double>(Nrac*7, 0.0));
vector<vector<double> > SurfaceAreaByLayer(Mlayers+1, vector<double>(Nrac*7, 0.0));
vector<vector<vector<double> > > OverByLayer(Mlayers+1, vector<vector<double> >(Nrac*7, vector<double>(Nrac*7, 0.0)));
for(int m=1; m<=10; m++) {
    char name[64];
    sprintf(name, "CoreStatsm%d.txt", m);
    fover = fopen(name, "r");
    if (fover != NULL) {
	int N;
	fscanf(fover, "%d", &N);
	if (N != Nrac*7) {
	    printf("Size mismatch %d %d in file '%s'\n", N, Nrac*7, name);
	    return 42;
	    }
	for(int y=0; y<N; y++)
	    if (fscanf(fover, "%lf %lf", &VolumeByLayer[m][y], &SurfaceAreaByLayer[m][y] ) != 2)
		return 42;
	for(int y=0; y<N; y++) {
	    for(int x=0; x<N; x++) {
		if (fscanf(fover, "%lf", &(OverByLayer[m][y][x]) ) != 1)
		    return 42;
		}
	    }
	fclose(fover);
	}
    }
// Compute the Peter's rule plots for each medulla layer
//for(int m=1; m<=10; m++) {
for(int m=1; m<0; m++) {  // Kill this for now
    char name[64];
    sprintf(name, "scatter%d/over", m);
    FILE *fp = fopen(name, "w");
    if (fp == NULL) {
	printf("Unable to open '%s'.\n", name);
	return 42;
	}
    for(int f=0; f<Nrac; f++) {
	int ff = f*7;
	for(int t=0; t<Nrac; t++) {
	    int tt = t*7;
	    double synapse[7];
	    double overlap[7];
	    for(int k=0; k<7; k++) {
		int from = subset_index[ff+k];
		int to   = subset_index[tt+k];
		where local(from, to, m);
		set<where>::iterator it = ss.find(local);
                int str = 0;
		if (it != ss.end())
		    str = it->n;
                synapse[k] = str;
                overlap[k] = OverByLayer[m][ff+k][tt+k]/10000.0;
                if (synapse[k] > 0 && overlap[k] < 1.0e-6)
		    printf("Odd in Peter calc - connection but no area: %d %d\n", ff+k, tt+k);
		}
	    //ProcessSynOverlaps(all_col_names[f], all_col_names[t], synapse, overlap, fp, m);
	    }
        }
    fclose(fp);
    }
// Autapse data
{
vector<int>    count(NumTypes, 0); // for averaging the percentages
vector<double>  sums(NumTypes, 0);
vector<int>    total_in (NumTypes, 0);
vector<int>    total_aut(NumTypes, 0);
for(int i=0; i<names.size(); i++) {
    if (!isdigit(names[i][0])) {
	int aut = from_to[i][i];
	int sum = 0;
	for(int k=0; k<names.size(); k++)
	    sum += from_to[k][i];
        int ty = TypeIndex(names[i]);
        double percent = double(aut)/sum*100.0;
	printf("Autapse (all) %s %d : %d of %d, %.2f%%\n", names[i], ty, aut, sum, percent);
        if (ty >= 0 /* && sum >= 50 */) {
	    count[ty]++;
            sums [ty] += percent;
            total_in [ty] += sum;
            total_aut[ty] += aut;
	    }
	}
    }
    const int BINS_PER_PERCENT = 20;
    vector<int> histo(3*BINS_PER_PERCENT,0);
    vector<int> bin_sum(3*BINS_PER_PERCENT, 0);
    vector<int> bin_aut(3*BINS_PER_PERCENT, 0);
    for(int k=0; k<NumTypes; k++) {
        double avg = sums[k]/max(1, count[k]); // average over all cell instances, or percentage
        double avg2= double(total_aut[k])/total_in[k] * 100.0;
	printf("Autapse (histo) %s %d %d %.2f %2f\n", CellTypes[k].name, k, count[k], avg, avg2);
        int bin = int(avg2*BINS_PER_PERCENT);
        if (bin >= 0 && bin < histo.size()) {
	    histo[bin]++;
            bin_sum[bin] += total_in [k];
            bin_aut[bin] += total_aut[k];
	    }
	}
    for(int k=0; k<histo.size(); k++) {
	double x1 = double(k)  / BINS_PER_PERCENT;
	double x2 = (k + 0.99) / BINS_PER_PERCENT;
        printf("%.3f %d\n%.3f %d\n", x1, histo[k], x2, histo[k]);
	}
    int cumul_sum = 0;
    int cumul_aut = 0;
    for(int k=0; k<histo.size(); k++) {
	cumul_sum += bin_sum[k];
	cumul_aut += bin_aut[k];
	printf("Cumulative including bin %d: %d of %d, %.3f%%\n", k, cumul_aut, cumul_sum, double(cumul_aut)/cumul_sum*100.0);
	}
}

FILE *f7 = fopen("f7.txt", "w");
if (f7 == NULL) {
    printf("Could not open 'f7.txt' for write.\n");
    exit(42);
    }
fover = fopen("scatter/over", "w");
if (fover == NULL) {
    printf("Could not open 'scatter/over' for write\n");
    exit(42);
    }
FILE *fz = fopen("scatter/overz", "w");
if (fz == NULL) {
    printf("Could not open 'scatter/overz' for write\n");
    exit(42);
    }
// print the individual tables.
printf("---Core Connectome, human readable--\n");
printf("The rows and columns (from left to right, and top to bottom) represent columns H,a,b,c,d,e,f\n");
printf("as described in the text.  The left table is the number of synapses that are pre-synaptic in the cell\n");
printf("of the row and post-synaptic to the cell in the appropriate column.\n");
printf("The right table is the contact area between the two cells, in square microns.\n");
printf("\n");
vector<vector<int> > total_over_cols(7, vector<int>(7, 0));
int naz = 0;       // number of all non-zero matrices
int npd = 0;       // number of pure diagonal, with all >0
vector<double> fitvar_x;
vector<double> fitvar_y;
vector<double> fitvar_e;   // for our estimate based on multipicity
vector<double> fitvar_m;   // multiplicity alone
for(int f = 0; f<Nrac; f++) {
    int ff = f*7;
    for(int t=0; t<Nrac; t++) {
	int tt = t*7;
        int count_nz = 0;
        int non_diag_nz = 0;
        double diag_sum  = 0.0;
	double diag_sum2 = 0.0;
        MeanStd mult_stats;
        for(int y=0; y<7; y++) {
	    for(int x=0; x<7; x++) {
		count_nz += int(ic[ff+y][tt+x] > 0);
		non_diag_nz += int(x != y && ic[ff+y][tt+x] > 0);
                total_over_cols[x][y] += ic[ff+y][tt+x];
		}
            diag_sum  += ic[ff+y][tt+y];
	    diag_sum2 += ic[ff+y][tt+y]*ic[ff+y][tt+y];
            // look up Pds/tbar ratio for diagonal connections
            int total_str = 0;
            int total_tbar= 0;
            where local(subset_index[ff+y], subset_index[tt+y], 0);
	    set<where>::iterator look = str_by_pair.find(local);
            if (look != str_by_pair.end())
		total_str = look->n;
            look = Tbar_by_pair.find(local);
            if (look != Tbar_by_pair.end())
		total_tbar = look->n;
            if (total_str > 0) {
                double rat = double(total_str)/total_tbar;
                printf("Variance of mult %d %d %.3f\n", total_str, total_tbar, rat);
                mult_stats.Element(rat);
		}
	    }
        if (mult_stats.HowMany() > 1)
	    printf("stats of mult - found %d, mean %.3f, std %.3f, var %.4f\n", mult_stats.HowMany(), mult_stats.Mean(), mult_stats.Std(), mult_stats.Var() );
        naz += int(count_nz == 0);
        npd += int(non_diag_nz == 0 && count_nz == 7);
	// print the rows and cols within the block.  Add a small message if all zero or purely diagonal.
	char num[32];
        sprintf(num, "(%d non-zero)", count_nz);
	printf("\nFrom %s to %s: %s\n", all_col_names[f], all_col_names[t], 
         (count_nz == 0) ? "(All zeros)" : (non_diag_nz == 0 ? "(Pure diagonal)": num) );
        for(int y=0; y<7; y++) {
	    for(int x=0; x<7; x++) 
		printf("%4d", ic[ff+y][tt+x]);
            printf("    ");
	    for(int x=0; x<7; x++) 
		printf("%8.2f", CoreOverlaps[ff+y][tt+x]/10000.0);  // convert to square microns
            printf("\n");
	    }
        if (f == t)
	    printf("Autapse %s: %d %.2f\n", all_col_names[f], int(diag_sum), diag_sum/7.0);
        if (diag_sum > 0) {
	    double diag_avg = diag_sum/7;
	    for(int y=0; y<7; y++) {
		if (diag_avg > 0)
		    printf("EQN: %.4f X%d = %d\n", diag_avg, y+1, ic[ff+y][tt+y] );
                }
	    double var = diag_sum2/7 - diag_avg*diag_avg;
	    double p = 0.5;  // guess at detection probability
            double approxN = diag_avg/p;
            double evar = approxN*p*(1-p); // expected variance of a binomial distribution, N terms, prob p;
            // look up ratio computed earlier by assuming missing connections really connected
            // we can do this over all connections of this type, or 7 colunm only
            double tbr    = 1.0;  // Tbar ratio
            double tbr_var = 0.0;  // and it's variance
            bool use_all = false;
            if (use_all) {
                string name = string(all_col_names[f]) + ":" + string(all_col_names[t]);
		map<string,double>::iterator lr =  from_to_ratio.find(name);
		if (lr == from_to_ratio.end())
		    printf("Bogons looking for %s\n", name.c_str() );
		else
		    tbr = lr->second / from_to_weight[name];
		}
	    else { // use just from 7 column, if available
		if (mult_stats.HowMany() >= 4) {
		    tbr = mult_stats.Mean();
                    tbr_var = mult_stats.Var();
		    }
	        }
            printf("%s - %s, Avg %.2f, Variance %.2f, measurement %.2f %s\n", 
	     all_col_names[f], all_col_names[t], diag_avg, var, evar, var < evar ? "Surprisingly good" : "");

            double mod = diag_avg*tbr*p + tbr_var*SQR(diag_avg/tbr) + diag_avg*(1-p);  // expected variance from our model
            printf("Mini: %.2f %.2f %.2f %.2f %.3f \"%s:%s\"\n", 
	    diag_avg, var, evar, mod, tbr, all_col_names[f], all_col_names[t]);
            if (diag_avg >= 1 && var >= diag_avg/3) {  // second term removes 3 outliers
                fitvar_x.push_back(log(diag_avg));
                fitvar_y.push_back(log(var));
                fitvar_e.push_back(log(mod));
                fitvar_m.push_back(log(tbr));
                }
            // now write into file designed for comparing columns
            for(int y=0; y<7; y++)
		fprintf(f7, "%4d,", ic[ff+y][tt+y]);
            double missing = ic[ff+1][tt+1] - diag_avg;
 	    fprintf(f7, "\"%s->%s\", %.2f\n", all_col_names[f], all_col_names[t], missing);
	    //
	    // Now create a small file for plotting scatter among the synpase/area ratio
	    double synapse[7];
	    double overlap[7];
            for(int y=0; y<7; y++) {
		synapse[y] = ic[ff+y][tt+y];
	        overlap[y] = CoreOverlaps[ff+y][tt+y] / 10000.0; // convert to square microns
		}
	    //ProcessSynOverlaps(all_col_names[f], all_col_names[t], synapse, overlap, fover, -1 /* all layers */);
	    }
        else {
            // all diagonal elements must be 0.  If the area is > 0, then
            // write to a special file, with a synpase/area of 0.1, so they
            // can appear on the bottom of the plot in a different color.
            MeanStd m;
	    for(int k=0; k<7; k++)
		m.Element(CoreOverlaps[ff+k][tt+k]/10000.0);
            if (m.Mean() > 0)
                fprintf(fz, "%.3f 0.1\n", m.Mean());
	    }
	}
    }
fclose(fover);
fclose(fz);
fclose(f7);
// Write a Matlab file for fits
FILE *fm = fopen("matlab.m", "w");
if (fm == NULL) {
    printf("Unable to open matlab file\n");
    return 42;
    }
{
vector<Point> temp(fitvar_x.size());
// here we sort the multiplicity data, then write it out
for(int i=0; i<fitvar_x.size(); i++)
    temp[i] = Point(fitvar_x[i],  fitvar_m[i]);
sort(temp.begin(), temp.end());
fprintf(fm, "x=[");
for(int i=0; i<fitvar_x.size(); i++)
    fprintf(fm, "%c%f", i==0? ' ' : ',', temp[i].x );
fprintf(fm, "]';\nm=[");
for(int i=0; i<fitvar_m.size(); i++)
    fprintf(fm, "%c%f", i==0? ' ' : ',', temp[i].y );
fprintf(fm, "]';\n");

// Now the same for the variance data
for(int i=0; i<fitvar_x.size(); i++)
    temp[i] = Point(fitvar_x[i],  fitvar_y[i]);
sort(temp.begin(), temp.end());
// don't need to re-write x
fprintf(fm, "y=[");
for(int i=0; i<fitvar_m.size(); i++)
    fprintf(fm, "%c%f", i==0? ' ' : ',', temp[i].y );
fprintf(fm, "]';\n");

// now for our estimate
for(int i=0; i<fitvar_x.size(); i++)
    temp[i] = Point(fitvar_x[i],  fitvar_e[i]);
sort(temp.begin(), temp.end());
// don't need to re-write x
fprintf(fm, "e=[");
for(int i=0; i<fitvar_m.size(); i++)
    fprintf(fm, "%c%f", i==0? ' ' : ',', temp[i].y );
fprintf(fm, "]';\n");
}
fclose(fm);

// Now do NR fits (this was the original
Fitab fvar(fitvar_x, fitvar_y);  // fvar.a is the constant, fvar.b the linear term.
printf("Fit to variance in log domain %f + %f*x\n", fvar.a, fvar.b);
const double BOT = 1.0;
const double TOP = 200.0;
double var_y1 = exp(fvar.a + fvar.b * log(BOT));
double var_y2 = exp(fvar.a + fvar.b * log(TOP));
printf("Fit to variance on log scale (%.1f,%f) to (%.1f, %f)\n", BOT, var_y1, TOP, var_y2);
// Now fit the estimate with multiplicity included
Fitab fvar_e(fitvar_x, fitvar_e);  // fvar.a is the constant, fvar.b the linear term.
var_y1 = exp(fvar_e.a + fvar_e.b * log(BOT));
var_y2 = exp(fvar_e.a + fvar_e.b * log(TOP));
printf("Est. with multiplicities to variance on log scale (%.1f,%f) to (%.1f, %f)\n", BOT, var_y1, TOP, var_y2);
// Now fit the multiplicity itself
Fitab fvar_m(fitvar_x, fitvar_m);  // fvar.a is the constant, fvar.b the linear term.
var_y1 = exp(fvar_m.a + fvar_m.b * log(BOT));
var_y2 = exp(fvar_m.a + fvar_m.b * log(TOP));
printf("fit to multiplicity on log scale (%.1f,%f) to (%.1f, %f)\n", BOT, var_y1, TOP, var_y2);
// Try direct proportion
vector<double> fitvar_s(fitvar_x.size(), 1.0);
Fitlin fvar2(fitvar_x, fitvar_y, fitvar_s, FitALinear);
fvar2.hold(1, 1.0);  // make the slope 1
fvar2.fit();
printf("Purely proportional %f %f\n", fvar2.a[0], fvar2.a[1]);
printf("Fit to variance on log scale (%.1f,%f) to (%.1f, %f)\n", BOT, BOT*exp(fvar2.a[0]), TOP, TOP*exp(fvar2.a[0]));

printf("Grand total by columns\n");
for(int f=0; f<7; f++) {
    for(int t=0; t<7; t++)
	printf(" %4d", total_over_cols[f][t]);
    printf("\n");
    }

// Do an unbiased search for correlations
vector<vector<double> > dd(Nrac*Nrac, vector<double>(7) );
MeanStd cols[7];
for(int f=0; f<Nrac; f++) {
    int ff = f*7;
    for(int t=0; t<Nrac; t++) {
	int tt = t*7;
        for(int k=0; k<7; k++) {
	    dd[f*Nrac+t][k] = ic[ff+k][tt+k];
	    cols[k].Element(ic[ff+k][tt+k]);
	    }
	}
    }
for(int k=0; k<7; k++) {
    printf("Sum of col %d = %6.1f, %6.2f%%\n", k, cols[k].Sum() , cols[k].Sum()/cols[0].Sum()*100.0);
    double rat = cols[0].Sum()/cols[k].Sum();
    for(int j=0; j<Nrac*Nrac; j++)
	dd[j][k] *= rat;
    }

SearchForCorrelation(dd, Nrac, all_col_names, '/', "NullRatio.txt");
SearchForCorrelation(dd, Nrac, all_col_names, '+', "NullAdd.txt");
printf("---End Core Connectome, human readable---\n");
// write the in-column component for each column
for(int k=0; k<7; k++) {
    printf("\nMatrix from <row> to <column> for column %c:\n", "Habcdef"[k]);
    printf("           ");
    for(int i=0; i<Nrac; i++)
	printf("%4d", i);
    printf("\n\n");
    for(int f=0; f<Nrac*7; f += 7) {
        printf("%2d %8s", f/7, all_col_names[f/7]);
	for(int t=0; t<Nrac*7; t += 7)
	    printf("%4d", ic[f+k][t+k]);
	printf("\n");
	}
    printf("\nMatrix to <row> from <column> for column %c:\n", "Habcdef"[k]);
    printf("           ");
    for(int i=0; i<Nrac; i++)
	printf("%4d", i);
    printf("\n\n");
    for(int t=0; t<Nrac*7; t += 7) {
        printf("%2d %8s", t/7, all_col_names[t/7]);
        for(int f=0; f<Nrac*7; f += 7)
	    printf("%4d", ic[f+k][t+k]);
	printf("\n");
	}
    }

printf("Grand:  %d matrices all zero, %d pure diagonal and all > 0\n", naz, npd);
// Write the core body IDs in json
FILE *fb = fopen("Core-bodies.json", "w");
if (fb == NULL) {
    printf("can't open core bodies for write.\n");
    return 42;
    }
fprintf(fb, "{\"data\": [\n");
for(int i=0; i<Nrac*7; i++)
    fprintf(fb, "%d%s", rev_map[subset_index[i]], i < Nrac*7-1 ? "," : "" );
fprintf(fb, "]}\n");
fclose(fb);

double grand_total_cons    = 0.0;  // sum of consistent
double grand_total_type1a  = 0.0;  // sum those with 1-3
double grand_total_type1b  = 0.0;  // sum those with 1-3
double grand_total_missing = 0.0;  // extrapolated missing amount
double grand_total_dc      = 0.0;  // sum of strengths used to compute missing; for double check.
vector<overlap> suspicious;        // has a pair of IDs + sorting, so use it to record suspicious connections (1a and 1b)
vector<double> histo_by_strength_consistent(50,0.0);
vector<double> histo_by_strength_incon_type1a(50,0.0);   // unexpected cell type
vector<double> histo_by_strength_incon_type1b(50,0.0);   // expected type, but extra connection
vector<double> histo_by_strength_incon_type2(50,0.0);
vector<int> histo_nzr(8,0); // histogram of this type to type count (outputs)
vector<int> histo_nzc(8,0); // histogram of      type to type inputs

// Data for table of errors by type
vector<double> out_err_type1a(Nrac, 0.0);
vector<double> out_err_type1b(Nrac, 0.0);
vector<double> out_err_type2 (Nrac, 0.0);
vector<double> out_err_type3 (Nrac, 0.0);  // type 3 is an autapse
vector<double> out_err_denom (Nrac, 0.0);  // denominator, total number of connections
vector<double>  in_err_type1a(Nrac, 0.0);
vector<double>  in_err_type1b(Nrac, 0.0);
vector<double>  in_err_type2 (Nrac, 0.0);
vector<double>  in_err_type3 (Nrac, 0.0);  // autapse
vector<double>  in_err_denom (Nrac, 0.0);  // denominator, total number of connections

int N_1a_cases = 0;
int N_1b_cases = 0;
int N_2_cases  = 0;

for(int f = 0; f<Nrac; f++) {  // f = 'from'
    printf("\n");
    int ff = f*7;
    for(int t=0; t<Nrac; t++) { // t = 'to'
	int tt = t*7;
	// sum the rows and cols within the block
	vector<int> row_sums(7,0);
	vector<int> col_sums(7,0);
        vector<int> row_nz(7,0);      // number of non-zero entries in each row
        vector<int> col_nz(7,0);      // number of non-zero entries in each icolumn
        double tot = 0.0;
        double trace = 0.0;
        for(int y=0; y<7; y++) {
	    for(int x=0; x<7; x++) {
		row_sums[y] += ic[ff+y][tt+x];
                col_sums[x] += ic[ff+y][tt+x];
                row_nz  [y] += int(ic[ff+y][tt+x] > 0);
                col_nz  [x] += int(ic[ff+y][tt+x] > 0);
                tot += ic[ff+y][tt+x];
                trace += (x == y) ? ic[ff+y][tt+x] : 0.0;  // trace of the matrix
		}
	    }
        out_err_denom[f] += tot;
         in_err_denom[t] += tot;
        bool nadr = true;   // No analysis done by row
        bool nadc = true;   // No analysis done by column
        // Look for type Ib errors.  First check where an output drives an extra cell 
        int HowMany;
        int mode = FindMode(row_nz, &HowMany);
        if ((mode == 1 && HowMany) >= 4 || (trace > 35 && trace/tot > 0.9) ) {  // If most of them connect to 1 thing, others are extra
            nadr = false;
	    N_1b_cases++;
	    for(int y=0; y<7; y++) {
		if (row_nz[y] > 1) {
                    // find the strength.  We want the sum of all but the largest
                    vector<int> sorted(7);
                    for(int k=0; k<7; k++)
			sorted[k] = ic[ff+y][tt+k];
                    sort(sorted.begin(), sorted.end());
		    int str = 0;
                    for(int k=0; k<6; k++)
			str += sorted[k];
                    const char *col = "Habcdef";
		    printf("Possible type 1b extra output between %s, col %c and %s, strength %d\n", all_col_names[f], col[y], all_col_names[t], str);
                    grand_total_type1b += str;
	            AddToStrengthHisto(histo_by_strength_incon_type1b, str);
                    out_err_type1b[f] += str;
		    // all the ones except the strongest are suspicious
		    for(int k=0; k<7; k++)
			if (ic[ff+y][tt+k] > 0 && ic[ff+y][tt+k] != sorted[6])
			    suspicious.push_back(overlap(subset_index[ff+y], subset_index[tt+k], 1.0));
		    }
		}
	    }
        // Look for type Ib errors.  Next check where a cell gets input from an extra cell 
        HowMany;
        mode = FindMode(col_nz, &HowMany);
        if ((mode == 1 && HowMany >= 4) || (trace > 35 && trace/tot > 0.9) ) {  // If most of them connect to 1 thing, others are extra
            nadc = false;
	    N_1b_cases++;
	    for(int x=0; x<7; x++) {
		if (col_nz[x] > 1) {
                    // find the strength.  We want the sum of all but the largest
                    vector<int> sorted(7);
                    for(int k=0; k<7; k++)
			sorted[k] = ic[ff+k][tt+x];
                    sort(sorted.begin(), sorted.end());
		    int str = 0;
                    for(int k=0; k<6; k++)
			str += sorted[k];
                    const char *col = "Habcdef";
		    printf("Possible type 1b extra  input between %s, col %c and %s, strength %d\n", all_col_names[f], col[x], all_col_names[t], str);
                    grand_total_type1b += str;
	            AddToStrengthHisto(histo_by_strength_incon_type1b, str);
                    in_err_type1b[t] += str;
		    // all the ones except the strongest are suspicious
		    for(int k=0; k<7; k++)
			if (ic[ff+k][tt+x] > 0 && ic[ff+k][tt+x] != sorted[6])
			    suspicious.push_back(overlap(subset_index[ff+k], subset_index[tt+x], 1.0));
		    }
		}
	    }
        // row_sum[y] = total strength of one source instance to all cells of target type
        // col_sum[x] = total strength of all inputs of source type to one instance of target
        // count how many non-zero entries
        int nzr = 0;  // non-zero rows.  Out of the 7 source cells, how many talk to at least one cell of the specified type
	int nzc = 0;  //                 Out of the 7 target cells, how many get  at    least one input from a cell of the specified type
        double rss = 0.0; // Row sum of sums
        double css = 0.0; // column sum of sums
        for(int i=0; i<7; i++) {
	    if (row_sums[i] > 0) nzr++;
            rss += row_sums[i];
	    if (col_sums[i] > 0) nzc++;
            css += col_sums[i];
	    }

	printf("From (out) %.8s to %.8s has %2d non-zeros (%.2f)  [", all_col_names[f], all_col_names[t], nzr, rss/max(nzr,1) );
	histo_nzr[nzr]++;
        double asnze = rss/nzr;  // average strength of non-zero entries
	if (nzr >= 1 && nzr <= 3) {
            nadr = false;
            N_1a_cases++;
	    AddToStrengthHisto(histo_by_strength_incon_type1a, asnze * nzr);
            out_err_type1a[f] += rss;
            grand_total_type1a += rss;
            for(int y=0; y<7; y++)
		for(int x=0; x<7; x++)
		    if (ic[ff+y][tt+x] > 0)
			suspicious.push_back(overlap(subset_index[ff+y], subset_index[tt+x], 0.0));
	    }
	if (nzr >= 4 && nzr <= 6) {
            nadr = false;
            N_2_cases++;
	    AddToStrengthHisto(histo_by_strength_incon_type2, asnze * (7-nzr));  // Amount missing
            out_err_type2[f] += asnze * (7-nzr);
	    grand_total_missing += asnze * (7-nzr);
            grand_total_dc      += asnze * nzr;
            // heuristic to find possible missing connections
            for(int k=0; k<7; k++) {
		if (ic[ff+k][tt+k] == 0) {
  		    LookForMissing(from_to, names, rev_map, subset_index[ff+k], subset_index[tt+k], asnze);
		    }
		}
	    }
        // Now compute autapses
        if (f == t) {
	    out_err_type3[t] += trace;
             in_err_type3[t]  += trace;
	    }
	if (nzr == 7) {
	    AddToStrengthHisto(histo_by_strength_consistent,   rss);
	    grand_total_cons += rss;
	    }
	//for(int i=0; i<7; i++) {
	    //printf("%4d", row_sums[i]);
	    //if (row_sums[i] > 0 && nzr >= 1 && nzr <= 3)
		//AddToStrengthHisto(histo_by_strength_incon_type1a, row_sums[i]);
	    //if (row_sums[i] > 0 && nzr >= 4 && nzr <= 6)
		//AddToStrengthHisto(histo_by_strength_incon_type2, row_sums[i]);  // actually, should set to missing
	    //if (row_sums[i] > 0 && nzr == 7)
		//AddToStrengthHisto(histo_by_strength_consistent,   row_sums[i]);
	    //}
	printf("]\n");

	printf("From ( in) %.8s to %.8s has %2d non-zeros (%.2f)  [", all_col_names[f], all_col_names[t], nzc, css/max(nzc,1) );
	histo_nzc[nzc]++;
        asnze = css/nzc;  // average strength of non-zero entries
	if (nzc >= 1 && nzc <= 3) {
            nadc = false;
            N_1a_cases++;
	    AddToStrengthHisto(histo_by_strength_incon_type1a, asnze * nzc);
            in_err_type1a[t] += css;
            grand_total_type1a += css;
            for(int y=0; y<7; y++)
		for(int x=0; x<7; x++)
		    if (ic[ff+y][tt+x] > 0)
			suspicious.push_back(overlap(subset_index[ff+y], subset_index[tt+x], 0.0));
	    }
	if (nzc >= 4 && nzc <= 6) {
            nadc = false;
	    N_2_cases++;
	    AddToStrengthHisto(histo_by_strength_incon_type2, asnze * (7-nzc));  // Amount missing
            in_err_type2[t] += asnze * (7-nzc);
	    grand_total_missing += asnze * (7-nzc);
            grand_total_dc      += asnze * nzc;
            // heuristic to find possible missing connections
            for(int k=0; k<7; k++) {
		if (ic[ff+k][tt+k] == 0) {
                    LookForMissing(from_to, names, rev_map, subset_index[ff+k], subset_index[tt+k], asnze);
		    }
		}
	    }
	if (nzc == 7) {
	    AddToStrengthHisto(histo_by_strength_consistent,   css);
            grand_total_cons += css;
	    }
	//for(int i=0; i<7; i++) {
	    //printf("%4d", col_sums[i]);
	    //if (col_sums[i] > 0 && nzc >= 1 && nzc <= 3)
		//AddToStrengthHisto(histo_by_strength_incon_type1a, col_sums[i]);
	    //if (col_sums[i] > 0 && nzc >= 4 && nzc <= 6)
		//AddToStrengthHisto(histo_by_strength_incon_type2, col_sums[i]);  // should be missing
	    //if (col_sums[i] > 0 && nzc == 7)
		//AddToStrengthHisto(histo_by_strength_consistent,   col_sums[i]);
	    //}
	printf("]\n");
        printf("Debugging cases:  %d type 1a,  %d type 1b, %d type 2\n", N_1a_cases, N_1b_cases, N_2_cases);
        if (tot > 0 && nadc && nadr)
	    printf("   ### No analysis done for this source-target pair\n");
	}
    }
double grand_sum = grand_total_type1a + grand_total_dc + grand_total_cons;
printf("Grand totals: %.2f extra 1a, %.2f extra 1b, %.2f missing 2, %.2f double check, %.2f consistent, %.2f sum\n",
 grand_total_type1a, grand_total_type1b, grand_total_missing, grand_total_dc, grand_total_cons, grand_sum);
printf("Grand percent type 1a %.2f%%, type 1b %.2f%%, type 2 %.2f%%\n", 
 grand_total_type1a/grand_sum*100.0, grand_total_type1b/grand_sum*100.0, grand_total_missing/grand_sum*100.0);
printf("Grand cases:  %d type 1a,  %d type 1b, %d type 2\n", N_1a_cases, N_1b_cases, N_2_cases);

FILE *ftab = fopen("tab.tex","w");
if (ftab == NULL) {
   printf("can't open tex file\n");
   return 42;
   }
fprintf(ftab, "\\documentclass[10pt]{article}\n"
              "\\usepackage{array}\n"
              "\\usepackage[margin=0.5in]{geometry}\n"
              "\\begin{document}\n");
fprintf(ftab, "\\begin{tabular}{ r r r r r r p{0.5cm} r r r r r p{0.5cm} r}\n");

// Write the error numbers by cell type.
vector<int> sums(10,0);
printf("---Table---\n");
printf("    Type  inputs    E+U   E+N     E-    Aut         outputs    E+U    E+N     E-    Aut\n");
fprintf(ftab, "Type &inputs&E+U & E+N  &  E-  & Aut  & &    outputs &  E+U &  E+N  &  E-  & Aut & & Total\\\\\n");
printf("---------------------------------------------------------------------------------------\n");
fprintf(ftab, "\\hline\n");
for(int i=0; i<Nrac; i++) {
    bool is_Mi1 = strcmp(all_col_names[i], "Mi1") == 0;  /* Mi1 does not count for autapse*/
    sums[0] += int(in_err_denom[i]); 
    sums[1] += int(in_err_type1a[i]);
    sums[2] += int(in_err_type1b[i]);
    sums[3] += int(in_err_type2 [i]);
    if (!is_Mi1)
	sums[4] += int(in_err_type3 [i]);
    sums[5] += int(out_err_denom[i]); 
    sums[6] += int(out_err_type1a[i]);
    sums[7] += int(out_err_type1b[i]);
    sums[8] += int(out_err_type2 [i]);
    if (!is_Mi1)
        sums[9] += int(out_err_type3 [i]);
    char star = (is_Mi1 ? '*' : ' ');
    printf("%8s %6d %6d %6d %6d %6d%c %14d %6d %6d %6d %6d%c %9d\n", all_col_names[i], 
     int(in_err_denom[i]), 
     int(in_err_type1a[i]),
     int(in_err_type1b[i]),
     int(in_err_type2 [i]),
     int(in_err_type3 [i]),
     star,
     int(out_err_denom[i]), 
     int(out_err_type1a[i]),
     int(out_err_type1b[i]),
     int(out_err_type2 [i]),
     int(out_err_type3 [i]),
     star,
     int(in_err_denom[i]) + int(out_err_denom[i]) );
    fprintf(ftab,"%8s&%6d&%6d&%6d&%6d&%6d%c& &%14d&%6d&%6d&%6d&%6d%c& & %9d\\\\\n", all_col_names[i], 
     int(in_err_denom[i]), 
     int(in_err_type1a[i]),
     int(in_err_type1b[i]),
     int(in_err_type2 [i]),
     int(in_err_type3 [i]),
     star,
     int(out_err_denom[i]), 
     int(out_err_type1a[i]),
     int(out_err_type1b[i]),
     int(out_err_type2 [i]),
     int(out_err_type3 [i]),
     star,
     int(in_err_denom[i]) + int(out_err_denom[i]) );
    }
printf("-------------------------------------------------------------------------\n");
fprintf(ftab, "\\hline\n");
printf("%8s %6d %6d %6d %6d %6d %15d %6d %6d %6d %6d\n", "Totals", 
 sums[0], sums[1], sums[2], sums[3], sums[4], sums[5], sums[6], sums[7], sums[8], sums[9]);
fprintf(ftab, "%8s&%6d&%6d&%6d&%6d&%6d& &%15d&%6d&%6d&%6d&%6d\\\\\n", "Totals", 
 sums[0], sums[1], sums[2], sums[3], sums[4], sums[5], sums[6], sums[7], sums[8], sums[9]);
printf("%8s       %6.2f%%%6.2f%%%6.2f%%%6.2f%%                %6.2f%%%6.2f%%%6.2f%%%6.2f%%\n", "Percent", 
 sums[1]/double(sums[0])*100.0, sums[2]/double(sums[0])*100.0, sums[3]/double(sums[0])*100.0, sums[4]/double(sums[0])*100.0,
 sums[6]/double(sums[5])*100.0, sums[7]/double(sums[5])*100.0, sums[8]/double(sums[5])*100.0, sums[9]/double(sums[5])*100.0);
fprintf(ftab, "%8s  & &  %6.2f\\%%&%6.2f\\%%&%6.2f\\%%&%6.2f\\%%  & & &%6.2f\\%%&%6.2f\\%%&%6.2f\\%%&%6.2f\\%%\\\\\n", "Percent", 
 sums[1]/double(sums[0])*100.0, sums[2]/double(sums[0])*100.0, sums[3]/double(sums[0])*100.0, sums[4]/double(sums[0])*100.0,
 sums[6]/double(sums[5])*100.0, sums[7]/double(sums[5])*100.0, sums[8]/double(sums[5])*100.0, sums[9]/double(sums[5])*100.0);
fprintf(ftab, "\\end{tabular}\n");
fprintf(ftab, "\\end{document}\n");
fclose(ftab);

// now percentage
printf("\n");

printf("---Table---\n");
printf("    Type  inputs    E+U    E+N     E-         outputs    E+U      E+N      E-\n");
printf("-----------------------------------------------------------------------------\n");
for(int i=0; i<Nrac; i++) {
    printf("%8s %6d %6.2f%% %5.2f%% %5.2f%% %15d %6.2f%% %6.2f%% %6.2f%%\n", all_col_names[i], 
     int(in_err_denom[i]), 
     in_err_type1a[i]/in_err_denom[i]*100.0,
     in_err_type1b[i]/in_err_denom[i]*100.0,
     in_err_type2 [i]/in_err_denom[i]*100.0,
     int(out_err_denom[i]), 
     out_err_type1a[i]/out_err_denom[i]*100.0,
     out_err_type1b[i]/out_err_denom[i]*100.0,
     out_err_type2 [i]/out_err_denom[i]*100.0);
    }
// Now write some data for looking for correlations.   Each line has cell volume, cell surface area, total I/O, then all the errors
printf("\nFor plotting\n");
for(int i=0; i<Nrac; i++) {
    double avg_vol = 0.0;
    double avg_area= 0.0;
    for(int k=i*7; k<i*7+7; k++) {
	avg_vol += CoreVolume[k];
	avg_area += CoreSurfaceArea[k];
	}
    avg_vol  /= 7*1.0e6; // average over columns; convert to cubic microns
    avg_area /= 7*1.0e4; // average over columns; convert to square microns
    printf("%.3f %.3f %6d %6d %6d %6d %6d %15d %6d %6d %6d\n", avg_vol, avg_area,
     int(in_err_denom[i]) + int(out_err_denom[i]), 
     int(in_err_denom[i]), 
     int(in_err_type1a[i]),
     int(in_err_type1b[i]),
     int(in_err_type2 [i]),
     int(out_err_denom[i]), 
     int(out_err_type1a[i]),
     int(out_err_type1b[i]),
     int(out_err_type2 [i]) );
     }

// Suspicious ids should be written out
FILE *fbook = fopen("Suspicious.json", "w");
if (fbook == NULL) {
    printf("could not open Suspicious.json for write\n");
    return 42;
    }
fprintf(fbook, "{\"data\": [\n");
bool wrote_any = false;
sort(suspicious.begin(), suspicious.end() );
for(int i=0; i<suspicious.size(); i++) {
    int n1 = suspicious[i].id1;
    int n2 = suspicious[i].id2;
    // if previous one exists, and it's the same, skip this one.
    if (i > 0 && n1 == suspicious[i-1].id1 && n2 == suspicious[i-1].id2)
	continue;
    int body_id_1 = rev_map[n1];
    int body_id_2 = rev_map[n2];
    if (i < 15)
        printf("Suspicious: %s (%d) -> %s (%d)\n", names[n1], body_id_1, names[n2], body_id_2);
    // OK, go find location of all of these
    int Num = 1;  // will number them
    for(int j=0; j<tbs.size(); j++) {
        if (tbs[j].body_id != body_id_1)
	    continue;
	for(int k=0; k<tbs[j].partners.size(); k++) {
	    if (tbs[j].partners[k].body_id != body_id_2)
		continue;
            int x = int(tbs[j].partners[k].pt.x);
            int y = int(tbs[j].partners[k].pt.y);
	    //printf("Found one! at %d %d %d\n", x, y, tbs[j].z);
            if (wrote_any)             // if there is a previous one
	        fprintf(fbook, ",\n"); // then we need a comma
            const char * stype = suspicious[i].area > 0 ? "Type 1b" : "Type 1a";
            fprintf(fbook, "{\"text\": \" SUSPICIOUS (%s) %s %s #%d\", \"body ID\": -1, \"location\": [%d, %d, %d]}", 
             stype, names[n1], names[n2], Num, x+2, y+2, tbs[j].partners[k].z); // offset by 2 pixels
            Num++;
            wrote_any = true;
	    }
	}
    }
fprintf(fbook, "\n],\n");
fprintf(fbook, "\"metadata\": {\n");
fprintf(fbook, "  \"username\": \"takemuras\",\n");
fprintf(fbook, "  \"software version\": \"1.9.6\",\n");
fprintf(fbook, "  \"description\": \"bookmarks\",\n");
fprintf(fbook, "  \"file version\": 1,\n");
fprintf(fbook, "  \"software revision\": \"4661M\",\n");
fprintf(fbook, "  \"computer\": \"takemuras-ws2.janelia.priv\",\n");
fprintf(fbook, "  \"date\": \"01-July-2014 14:17\",\n");
fprintf(fbook, "  \"session path\": \"/groups/flyem/data/MLOLP-Z1112-19R/align1/layers/session/shinomiyak.AxonTrace1.20140602\",\n");
fprintf(fbook, "  \"software\": \"Raveler\"\n");
fprintf(fbook, "  }\n");
fprintf(fbook, "}\n");

fclose(fbook);

printf("Histograms for outs:");
for(int i=0; i<=7; i++)
    printf("%4d", histo_nzr[i]);
printf("\n");
printf("Histograms for ins:");
for(int i=0; i<=7; i++)
    printf("%4d", histo_nzc[i]);
printf("\n");

fp = fopen("hc", "w");
if (fp == NULL) {
    printf("Cannot write 'hc'\n");
    return 42;
    }
for(int i=0; i<histo_by_strength_consistent.size(); i++) {
    double bin_min = pow(1.2,i)/7;
    double bin_max = bin_min*1.2;
    // create a flat top bin
    fprintf(fp,"%f %.3f\n", bin_min+0.001, histo_by_strength_consistent[i]);
    fprintf(fp,"%f %.3f\n", bin_max-0.001, histo_by_strength_consistent[i]);
    }
fclose(fp);
fp = fopen("h1a", "w");
if (fp == NULL) {
    printf("Cannot write 'h1a'\n");
    return 42;
    }
const double offset = 1.02; // for offset on plot
for(int i=0; i<histo_by_strength_consistent.size(); i++) {
    double bin_min = pow(1.2,i)/7 / offset;
    double bin_max = bin_min*1.2;
    // create a flat top bin
    fprintf(fp,"%f %.3f\n", bin_min+0.001, histo_by_strength_incon_type1a[i]);
    fprintf(fp,"%f %.3f\n", bin_max-0.001, histo_by_strength_incon_type1a[i]);
    }
fclose(fp);
fp = fopen("h1b", "w");
if (fp == NULL) {
    printf("Cannot write 'h1b'\n");
    return 42;
    }
for(int i=0; i<histo_by_strength_consistent.size(); i++) {
    double bin_min = pow(1.2,i)/7 / offset / offset;
    double bin_max = bin_min*1.2;
    // create a flat top bin
    fprintf(fp,"%f %.3f\n", bin_min+0.001, histo_by_strength_incon_type1b[i]);
    fprintf(fp,"%f %.3f\n", bin_max-0.001, histo_by_strength_incon_type1b[i]);
    }
fclose(fp);
fp = fopen("h2", "w");
if (fp == NULL) {
    printf("Cannot write 'h2'\n");
    return 42;
    }
for(int i=0; i<histo_by_strength_consistent.size(); i++) {
    double bin_min = pow(1.2,i)/7 * offset;
    double bin_max = bin_min*1.2;
    // create a flat top bin
    fprintf(fp,"%f %.3f\n", bin_min+0.001, histo_by_strength_incon_type2[i]);
    fprintf(fp,"%f %.3f\n", bin_max-0.001, histo_by_strength_incon_type2[i]);
    }
fclose(fp);

vector<int> WhatType(M);
for(int i=0; i<M; i++)
    WhatType[i] = TypeIndex(names[i]);
// will be vector of strongest connections seen
vector<vector<int> > Strongest(NumTypes, vector<int>(NumTypes,0));
vector<vector<int> > OutProfile(M, vector<int>(NumTypes, 0));
vector<vector<int> > InProfile(M, vector<int>(NumTypes, 0));

// Create lists of unnamed cells with most connections.  One table will be by biggest single connection, other by sum of all connections.
//
vector<sort_down> from_list;  // where the from cell is unnamed
vector<sort_down>   to_list;  // where the 'to' cell is unnamed
vector<sort_down> both_list;
vector<sort_down> weak_named; // weak connections between named cell
const int WEAK = 4;
for(int i=0; i<M; i++) {
    bool i_is_named = !isdigit(names[i][0]);
    for(int j=0; j<M; j++) {
	if (from_to[i][j] > 0) {
            int strength = from_to[i][j];
            bool j_is_named = !isdigit(names[j][0]);
	    sort_down s(i, j, double(strength));
            if (!i_is_named)  // the from cell is unnamed
	        from_list.push_back(s);
            if (!j_is_named)  // the 'to' cell is unnamed
                to_list.push_back(s);
            if ((!i_is_named) && !j_is_named)
                both_list.push_back(s);
            if (i_is_named && j_is_named && strength <= WEAK)
		weak_named.push_back(s);
            int ti = WhatType[i];
            int tj = WhatType[j];
            Strongest[ti][tj] = max(Strongest[ti][tj], int(strength));
            OutProfile[i][tj] = max(OutProfile[i][tj], strength);  // strongest connection from cell i to a cell of type j
            InProfile[j][ti]  = max( InProfile[j][ti], strength);
            }
        }
    }
printf("%d connections from unnamed cells\n", from_list.size());
sort(from_list.begin(), from_list.end());
for(int i=0; i<100 && i < from_list.size() && from_list[i].val >= 10; i++) { // print the biggest
    printf("# syn %5d: %s -> %s\n", int(from_list[i].val), names[from_list[i].index], names[from_list[i].index2] );
    }
// Now for connections to unnamed cells
printf("%d connections to unnamed cells\n", to_list.size());
sort(to_list.begin(), to_list.end());
for(int i=0; i<100 && i < to_list.size() && to_list[i].val >= 10; i++) { // print the biggest
    printf("# syn %5d: %s -> %s\n", int(to_list[i].val), names[to_list[i].index], names[to_list[i].index2] );
    }
printf("%d between unnamed cells\n", both_list.size());
sort(both_list.begin(), both_list.end());
for(int i=0; i<100 && i < both_list.size() && both_list[i].val >= 10; i++) { // print the biggest
    printf("# syn %5d: %s -> %s\n", int(both_list[i].val), names[both_list[i].index], names[both_list[i].index2] );
    }
printf("Total of %d weak (s <= %d) connections between named cells.\n", weak_named.size(), WEAK);

fbook = fopen("WeakBook.json", "w");
if (fbook == NULL) {
    printf("could not open WeakBook.json for write\n");
    return 42;
    }
fprintf(fbook, "{\"data\": [\n");
wrote_any = false;
for(int i=0; i<weak_named.size(); i++) {
    int n1 = weak_named[i].index;
    int n2 = weak_named[i].index2;
    int body_id_1 = rev_map[n1];
    int body_id_2 = rev_map[n2];
    if (i < 15)
        printf("%s (%d) -> %s (%d)\n", names[n1], body_id_1, names[n2], body_id_2);
    // OK, go find location of all of these
    int Num = 1;  // will number them
    for(int j=0; j<tbs.size(); j++) {
        if (tbs[j].body_id != body_id_1)
	    continue;
	for(int k=0; k<tbs[j].partners.size(); k++) {
	    if (tbs[j].partners[k].body_id != body_id_2)
		continue;
            int x = int(tbs[j].partners[k].pt.x);
            int y = int(tbs[j].partners[k].pt.y);
	    //printf("Found one! at %d %d %d\n", x, y, tbs[j].z);
            if (wrote_any)             // if there is a previous one
	        fprintf(fbook, ",\n"); // then we need a comma
            fprintf(fbook, "{\"text\": \" WEAK %s %s #%d\", \"body ID\": -1, \"location\": [%d, %d, %d]}", 
             names[n1], names[n2], Num, x+2, y+2, tbs[j].partners[k].z); // offset by 5 pixels
            Num++;
            wrote_any = true;
	    }
	}
    }
fprintf(fbook, "\n],\n");
fprintf(fbook, "\"metadata\": {\n");
fprintf(fbook, "  \"username\": \"takemuras\",\n");
fprintf(fbook, "  \"software version\": \"1.9.6\",\n");
fprintf(fbook, "  \"description\": \"bookmarks\",\n");
fprintf(fbook, "  \"file version\": 1,\n");
fprintf(fbook, "  \"software revision\": \"4661M\",\n");
fprintf(fbook, "  \"computer\": \"takemuras-ws2.janelia.priv\",\n");
fprintf(fbook, "  \"date\": \"01-July-2014 14:17\",\n");
fprintf(fbook, "  \"session path\": \"/groups/flyem/data/MLOLP-Z1112-19R/align1/layers/session/shinomiyak.AxonTrace1.20140602\",\n");
fprintf(fbook, "  \"software\": \"Raveler\"\n");
fprintf(fbook, "  }\n");
fprintf(fbook, "}\n");

fclose(fbook);

// Now, for every unnamed cell, see what types we might be able to rule out.
// We will do this in the order of total number of connections, to get the most important ones first.
vector<sort_down> ByConnNumber;
for(int i=0; i<M; i++) {
    if (!isdigit(names[i][0]))  // already named?
	continue;
    int sum = 0;
    for(int k=0; k<NumTypes-1; k++) // don't use "Unknown"
	sum += (InProfile[i][k] + OutProfile[i][k]);
    ByConnNumber.push_back(sort_down(i, double(sum)));
    }
sort(ByConnNumber.begin(), ByConnNumber.end() );
for(int c=0; c<ByConnNumber.size(); c++) {
    int i = ByConnNumber[c].index;
    vector<sort_down> CouldBe;
    int Nlikely = 0;
    for(int k=0; k < NumTypes-1; k++) {
	double odds = CouldBeOdds(k, InProfile[i], OutProfile[i], Strongest);
        if (odds > 0.5)
	    Nlikely++;
        CouldBe.push_back(sort_down(k, odds));
	}
    if (Nlikely > 3)   // If too many good matches, skip (for now)
	continue;
    sort(CouldBe.begin(), CouldBe.end() );
    // print top 3 matches with their score.. If we have exactly one perfect match, mark with stars
    bool one_perfect = CouldBe[0].val > 0.999 && CouldBe[1].val < 0.999;
    const char *prefix =  one_perfect ? "MAYBE:" : "";
    printf("%s\n%sCell %s could be: ", prefix, prefix, names[i]);
    for(int j=0; j<3; j++) {
	int k = CouldBe[j].index;
	printf("(%s %f)", CellTypes[k].name, CouldBe[j].val);
        }
    printf("\n");
    // print the evidence
    printf("%s  Outs:", prefix);
    for(int k=0; k<NumTypes-1; k++) {
	if (OutProfile[i][k] > 0)
	    printf(" (%d, %s)", OutProfile[i][k], CellTypes[k].name);
	}
    printf("\n");
    printf("%s  Ins:", prefix);
    for(int k=0; k<NumTypes-1; k++) {
	if (InProfile[i][k] > 0)
	    printf(" (%d, %s)", InProfile[i][k], CellTypes[k].name);
	}
    printf("\n");
    }
printf("All done!\n");
return 0;
}
