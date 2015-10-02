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


FILE * logfile;

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


vector<NeuroTr> NeuroTransmitters;

// First column is cell type, second is color, third is index in Single-column paper, if defined.
CellType CellTypes[] = {
 CellType("PPL1-05",   87, true, 17),
 CellType("PPL1-06",   85, true, 17),
 CellType("MBON-06",  95, true,  5),
 CellType("MBON-07",  97, true,  5),
 CellType("MBON-11",  99, true,  5),
 CellType("MBON-14",  93, true,  5),
 CellType("MBON-18",  91, true,  5),
 CellType("MBON-19",  89, true,  5),
 CellType("MBON-X",   90, true,  5),
 CellType("KC-prime",  70    ),  // Added
 CellType("KC-c",   71, true, 21),
 CellType("KC-p",   72, true,  9),
 CellType("KC-s",   73, true,  9),
 CellType("KC-a",   74, true, 11),
 CellType("Other",   55, true,  6),   // Added
 CellType("Irr",     6, true, 13),
 CellType("alpha",   24, true    ),   // columnar, but not in old
 CellType("non",  11, true,  2),
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
 vector<MeanStd> &stats, 
 vector<int> &Tbar_count)
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
    fprintf(logfile, "Expansion is %d to %d\n", rows[i].size(), s.size() );
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
double margin =  1.0;  // use a bigger value for Mi1 to avoid end of table
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
    fprintf(logfile, "Column %d, %d are used, smallest %.2f\n", i, num, smallz);
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


//-------------------------------------- Now write the table of ALL connections, not just those in table
total = 0;                          // Now re-do total to include all
for(int i=0; i<NumTypes; i++) {
    if (stats[i].HowMany() > 0) {
        double mean, std;
        stats[i].Stats(mean, std);
        total += int(stats[i].HowMany() * mean + 0.5);  // should be int, but round...
        fprintf(logfile, " Stats %d, %f\n", stats[i].HowMany(), mean);
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
    fprintf(logfile, "sa[%d] %s %s %d %f\n", i, bname, suffix, sa[i].index, sa[i].val);
sort(sa.begin(), sa.end()-1);
for(int i=0; i<NumTypes; i++)
    fprintf(logfile, "sa[%d] %s %s %d %f\n", i, bname, suffix, sa[i].index, sa[i].val);
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

fprintf(logfile, "RF: ------------------ %s ------------------\n", names[row_ids[0]] );
// Using the previous table, write a set of commands for computing receptive fields
for(int i=0; i<row_ids.size(); i++) {
    fprintf(logfile, "AVG %s", NoSpace(names[row_ids[i]]));
    // Now for the 5 most numerous input types, print contribution
    for(int k=0; k<5; k++) {
	int ty_look_for = any[k].index;
        int sum = 0;
        for(int j=0; j < rows[i].size(); j++) {
	    int ty = TypeIndex(names[rows[i][j]]);
            if (ty == ty_look_for)
		sum += dat[i][j];
	    }
	fprintf(logfile, " %d %s-from-%s", sum, NoSpace(names[row_ids[i]]), CellTypes[ty_look_for].name);
	}
    fprintf(logfile, "\n");
    // now figure out how each of them got here
    for(int k=0; k<5; k++) {
	int ty_look_for = any[k].index;
        fprintf(logfile, " AVG %s-from-%s", NoSpace(names[row_ids[i]]), CellTypes[ty_look_for].name);
        for(int j=0; j < rows[i].size(); j++) {
	    int ty = TypeIndex(names[rows[i][j]]);
            if (ty == ty_look_for){
		fprintf(logfile, " %d %s", dat[i][j], NoSpace(names[rows[i][j]]) );
            }
	    }
	fprintf(logfile, "\n");
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

// breaks a list up into comma separated parts.  Allocates results so the argument does not need to be saved.
vector<const char *> CommaList(char *list) {
vector<const char *> rslt;
for(char *p = strtok(list, ","); p != NULL; p = strtok(NULL, ",") )
    rslt.push_back(strdup(p));
return rslt;
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
	    fprintf(logfile, "Mike: %6d from %6d at x,y=%.1f %.1f %s\n", what, tbs[i].body_id, pt.x, pt.y,
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

fprintf(logfile, "Sorted by source id\n");
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
        fprintf(logfile, "From %6d at %.1f %.1f\n", p[j].body_id, p[j].pt.x, p[j].pt.y);
	sumx += p[j].pt.x;
        sumy += p[j].pt.y;
	}
    fprintf(logfile, "\n");
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


// Trace back inputs starting with Body ID 'root'.
// Connectivity matrix is 'from_to'.
// 'wanted' is the cell types you want to trace.
// 'trace' is false for those cells you do not want to trace back from (typically inputs)
void CircuitBacktrace(
 int root, vector<vector<unsigned short> > &from_to, vector<char *> &names, vector<bool> wanted, vector<bool> trace, int threshold)
{
fprintf(logfile, "Michael:\n");
vector<bool> cell_used(names.size(), false);  // to avoid loops
vector<int> look_at;
look_at.push_back(root);  // for now, start with one T4
vector<int> more;

for(int hops=0; look_at.size() > 0; look_at = more) {
    fprintf(logfile, "----------- Cells %d hops from target----------\n", hops++);
    more.clear();
    vector<bool> req(names.size(), false);  // requested already on this pass?
    for(int k=0; k<look_at.size(); k++) 
	req[look_at[k]] = true;             // already being considered
    for(int k=0; k<look_at.size(); k++) {
        int n = look_at[k];
        fprintf(logfile, " Cell: %s\n", names[n] );
        cell_used[n] = true;
        // add all inputs from cells of wanted types.
        for(int i=0; i<names.size(); i++) {
            int str = from_to[i][n];
            if (str <= threshold)
		continue;
            int ti = TypeIndex(names[i]);
	    if(wanted[ti]) {
		fprintf(logfile, "  Input from %-12s, %3d synapses\n", names[i], str );
		// if we want to trace, and not requested already, add to list    
		if (trace[ti] && (!cell_used[i]) && !req[i]) {
		    more.push_back(i);
		    req[i] = true;  // now been requested
		    }
		}
	    }
	}
    }
fprintf(logfile, "------------------------ End of Circuit --------------------------\n");
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
fprintf(logfile, "Opened '%s', got %d entries\n", name, N);
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
 vector<vector<unsigned short> > &from_to, map<int, int> &rev_map) {
vector<int> outs;
for(int i=0; i<names.size(); i++)
    if (strncmp(names[i],name,strlen(name)) == 0) {
	fprintf(logfile, "Neuron %s counts for %d and %.4f\n", names[i], (1 << outs.size()), pow(0.01, outs.size()+1) );
	outs.push_back(i);
	}
fprintf(logfile, "found %d of %s\n", outs.size(), name );
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
	fprintf(logfile, "\n");
    fprintf(logfile, "Shinya: %4d %9d %12s %.4f\n", i, rev_map[KCs[sa[i].index]], names[KCs[sa[i].index]], sa[i].val );
    }
}

int main(int argc, char **argv)
{

char *logfilename = strdup("logfile.txt");
if (argc > 1) 
{
    logfilename = argv[1];
}
logfile = fopen(logfilename, "w");

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
fprintf(logfile, "Starting to read file %s\n", fname); fflush(stdout);
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
fprintf(logfile, "JSON body file read, starting to parse.\n"); fflush(stdout);
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
fprintf(logfile, "Starting to read file %s\n", fname); fflush(stdout);
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
fprintf(logfile, "Synapse annotation file read, starting to parse.\n"); fflush(stdout);
parsingSuccessful = reader.parse(&huge[0], syn_root );
if ( !parsingSuccessful )
{
    // report to the user the failure and their locations in the document.
    std::cout  << "Failed to parse configuration\n"
               << reader.getFormatedErrorMessages();
    exit(42);
    }

fprintf(logfile, "parsed file...\n");
Json::Value md = syn_root["metadata"];
std::string des = md["description"].asString();
int         ver = md["version"].asInt();

fprintf(logfile, "Description: '%s', version %d\n", des.c_str(), ver );

// OK, read and parsed both files.  Start by going through BODY file.  Need to consolidate all CT1
// neurons.  This is a special case where they look like separate body IDs, but we know (from external
// evidence) they are all the same neuron, one that joins outside our valume.  Create a set of ints,
// CT1s, of all the CT1 body IDS.

set<int> CT1s;
bool combine_CT1s = false;
Json::Value data = bdy_root["data"];
if (combine_CT1s) {
for(unsigned int i=0; i<data.size(); i++) {
    //fprintf(logfile, "getting entry %d\n", i);
    Json::Value v = data[i];
    // read the different things that can be in 'data'.  For now only "body ID" and "name"
    int t = v["body ID"].asInt();
    fprintf(logfile, "Body ID %d\n", t);
    Json::Value part = v["name"];
    if (!part.isNull()) {
	string name = part.asString();
	fprintf(logfile, "Name %s\n", name.c_str() );
        if (strncmp(name.c_str(), "CT1", 3) == 0) {
	    fprintf(logfile, "Found %s, ID %d\n", name.c_str(), t);
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
    //fprintf(logfile, "getting entry %d\n", i);
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
fprintf(logfile, "Found %d unique IDs, %d T-bars, %d PSDs\n", name_map.size(), Tcount, Scount);

//---------------------- now go through the JSON body file again
vector<char> leaves(name_map.size(), 'U');  // 'U' for unknown, 't' for tangential, 'o' for orphan
fd = fopen("AllNamedBodies.json","w");  // write a file of scripts to get volumes and areas
if (fd == NULL) {
    printf("Could not open script file 'AllNamedBodies.json'\n");
    return 42;
    }
fprintf(fd, "{\"data\":[\n");
int NPrint = 0;
int Ncount = 0;
data = bdy_root["data"];
for(unsigned int i=0; i<data.size(); i++) {
    //printf("getting entry %d\n", i);
    Json::Value v = data[i];
    // read the different things that can be in 'data'.  For now only "body ID" and "name"
    int t = v["body ID"].asInt();
    fprintf(logfile, "Body ID %d\n", t);
    // Check to make sure it does not exist, then add it to the map
    bool in_the_map = name_map.find(t) != name_map.end();
    //if (!in_the_map)
	//printf("ID %d has no synapses??\n", t);
    Json::Value part = v["name"];
    if (!part.isNull()) {
	string name = part.asString();
        fprintf(logfile, "Body name %s\n", name.c_str() );
        if (!in_the_map) {
	    fprintf(logfile, "Named body %s has no synapses??\n", name.c_str() );
	    continue;
	    }
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
fprintf(logfile, "Found %d unique IDs, %d T-bars, %d PSDs, %d named bodies\n", name_map.size(), Tcount, Scount, Ncount);
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
    fprintf(logfile, "Reading %d sizes and areas\n", N);
    for (int c = getc(fd); c!=EOF; c = getc(fd)) {
	if (isdigit(char(c))) {
	    ungetc(char(c), fd);
            int id;
            double vol, area;
            fscanf(fd, "%d %lf %lf", &id, &vol, &area);
	    fprintf(logfile, "ID %d, vol %.2f, area %.2f\n", id, vol/1e6, area/1e4);
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
        fprintf(logfile, "Got '%s'\n", buffer); fflush(stdout);
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
vector< vector<unsigned short> > from_to(M, vector<unsigned short>(M,0));
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
	fprintf(logfile, "TB convergent %d, flagged %d, roi %d\n", tb.convergent, tb.flagged, tb.roi);
        if (CT1s.find(tb.body_id) != CT1s.end())   // If any of the CT1s, replace by first
	    tb.body_id = *(CT1s.begin());
        if (name_map.find(tb.body_id) == name_map.end()) {
	    fprintf(logfile, "No %d for TBar in map\n", tb.body_id);
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
            fprintf(logfile, "  Flag for psd is %d\n", prt.flagged);
	    if (CT1s.find(prt.body_id) != CT1s.end())   // If any of the CT1s, replace by first
		prt.body_id = *(CT1s.begin());
            if (name_map.find(prt.body_id) == name_map.end()) {
	        fprintf(logfile, "No %d for partner in map\n", prt.body_id);
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
                    fprintf(logfile, "Autapse:  %s Tbar at %d %d %d, psd at %d %d %d\n", name, 
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
        fprintf(logfile, " --- Name %s, pre %d/%d,   post %d/%d\n", names[i], pre, num_pre, post, num_post);
        if (isdigit(names[i][0]) &&
         70 <=  pre &&   pre <= 130 && 40 <=  num_pre &&   num_pre <= 80 &&
	 130 <= post && post <= 250 && 70 <= num_post && num_post <= 130) { //looks like a KC
	    char new_name[128];
            if (num_post >= 100)
	        sprintf(new_name, "KCs-guess-%d", frag_id++);
            else
	        sprintf(new_name, "KCc-guess-%d", frag_id++);
	    fprintf(logfile, "Cell with name %s had %d mentions and is renamed to %s\n", names[i], mentioned[i], new_name);
	    names[i] = strdup(new_name);  // memory leak, but who cares
            }
	}
    }
// If a cell is mentioned more than MENT_THRESHOLD times, but has no name, make it "Unknown".
for(int i=0; i<names.size(); i++) {
    if (mentioned[i] >= MENT_THRESHOLD && isdigit(names[i][0])) {
	char new_name[128];
        sprintf(new_name, "Fragment-%d", frag_id++);
	fprintf(logfile, "Cell with name %s had %d mentions and is renamed to %s\n", names[i], mentioned[i], new_name);
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
	fprintf(logfile, "GUESS c=%d s=%d p=%d\n",c, s, p);
        char *old = names[i];
	if (c > 10*(s+p)) {names[i] = strdup("KC-c"); changes++;}
	if (s > 10*(c+p)) {names[i] = strdup("KC-s"); changes++;}
	if (p > 10*(c+s)) {names[i] = strdup("KC-p"); changes++;}
        if (names[i] != old)
	    fprintf(logfile, "Changed BodyId %d to %s\n", rev_map[i], names[i]);
	}
    fprintf(logfile, "---GUESS: %d changes\n", changes);
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
    fprintf(logfile, "PART: There are %d Kenyon Cells\n", KCs.size() );
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
	fprintf(logfile, "Component starting at %d\n", i);
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
        fprintf(logfile, "%6d %8d %12s ", i, rev_map[KCs[i]], names[KCs[i]]);
        vector<char> cs;
	for(int k=0; k<rslts.size(); k++) {
	    char c = FindMembership(i, rslts[k], KCs, names);
	    fprintf(logfile, "%c", c);
            cs.push_back(c);
	    }
	fprintf(logfile, "\n");
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
		fprintf(logfile, "--------Susicious classification.  Body ID %d is labeled %s but best matches are ", rev_map[KCs[i]], name);
                for(int k=0; k<cs.size(); k++)
		    fprintf(logfile, "%c",cs[k]);
		fprintf(logfile, "\n");
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
		fprintf(logfile, "--------It's unanimous.  Body ID %d is labeled %s but looks like KC-%c\n", rev_map[KCs[i]], name, cl);
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
    fprintf(logfile, "\nGroup %d, adding %d\n", group_id, i);
    vector<int> gp;
    while (!st.empty()) {
	int curr = st.top();
	st.pop();
	fprintf(logfile, " Just popped %d\n", curr);
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
			fprintf(logfile, "Tbar %d, [%d] <- %.2f -> Tbar %d, [%d] \n", i, j, d, k, m);
			tbs[k].gid = group_id;
			st.push(k);
			fprintf(logfile, " Pushing %d\n", k);
		        break;
			}
		    }
		}
	    }
        if (st.size() == old)
	    fprintf(logfile, "None found within range.  Best = %.2f\n", best);
	}
    groups.push_back(gp);
    group_id++;
    }
// Print them out.
for(int i=0; i<groups.size(); i++) {
    fprintf(logfile, "\nGroup %d\n", i);
    map<int,int> m;
    for(int j=0; j<groups[i].size(); j++) {
        int t = groups[i][j];
	fprintf(logfile, "  Tbar id %8d%c%c -> ", tbs[t].body_id, tbs[t].convergent ? '*': ' ', tbs[t].flagged ? '+' : ' ');
	for(int k=0; k < tbs[t].partners.size(); k++) {
            int n = tbs[t].partners[k].body_id;
	    fprintf(logfile, "%8d%c ", n, tbs[t].partners[k].flagged ? '+' : ' ');
            if (m.find(n) == m.end() )
		m[n] = 1;
	    else
		m[n] = m[n]+1;
	    }
	fprintf(logfile, "\t(%5d, %5d, %5d)", int(tbs[t].pt.x), int(tbs[t].pt.y), tbs[t].z);
	fprintf(logfile, "\n");
	}
    map<int,int>::iterator mi;
    for(mi = m.begin(); mi != m.end(); mi++)
	fprintf(logfile, "%9d %2d%c\n", mi->first, mi->second, mi->second == groups[i].size() ? '*' : ' ');
   }
// Zipursky code here ------------------------------------------------------------
int tiL2 = TypeIndex("L2");
fprintf(logfile, "Zipurksy, looking for %d\n", tiL2);
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
	    fprintf(logfile, "No %d for TBar in map\n", tb.body_id);
	    continue;
	    }
	char *name = names[name_map[tb.body_id]];
        if (TypeIndex(name) != tiL2)
	    continue;
	tb.pt.x = lo[0u].asDouble();
	tb.pt.y = lo[1u].asDouble();
        tb.z    = lo[2u].asInt();
	fprintf(logfile, " Tb-bar on %s, partners [", name );
        vector<string> pnames;
	for(int j=0; j<part.size(); j++) {
	    Partner prt;
	    prt.confidence = part[j]["confidence"].asDouble();
	    prt.body_id    = part[j]["body ID"].asInt();
            if (name_map.find(prt.body_id) == name_map.end()) {
	        fprintf(logfile, "No %d for partner in map\n", prt.body_id);
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
	    fprintf(logfile, "%12s", pnames[j].c_str());
            if (j > 0 && pnames[j].compare(pnames[j-1]) == 0)
		ndup++;
	    }
	fprintf(logfile, " ] %d dup\n", ndup);
	}
    }
fprintf(logfile, "End Zipursky\n");
// End Zipursky ------------------------------------------------------------------
// Percentage complete code
fprintf(logfile, "Completeness by PSDs identified\n");
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
    fprintf(logfile, " Cell %15s: %4d known, %4d unknown, %.2f%% complete\n", names[from], known, unknown, known/double(known+unknown)*100.0);
    }
fprintf(logfile, "End of Completeness\n");
// End of completeness -----------------------------------------------------------
bool slow = true;
if (slow) {
fprintf(logfile, "Unknown cells:\n");
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
    fprintf(logfile, " Cell %15s: %4d known in, %4d unknown in, %4d known out, %4d unknown out, %4d total\n", names[i], ink, inu, outk, outu, ink+inu+outk+outu);
    fflush(stdout);
    sort_down s(i, ink + outk);
    sd.push_back(s);
    }
fprintf(logfile, "End of Unknown.  %d candidates\n", sd.size() );
fflush(stdout);
sort(sd.begin(), sd.end());
for(int j=0; j<25 && j < sd.size(); j++) {
    int i = sd[j].index;
    fprintf(logfile, "Cell %s has %d known connections\n", names[i], int(sd[j].val) );
    vector<int> ins_by_type(NumTypes,0);
    vector<int> out_by_type(NumTypes,0);
    vector<int> str_ins_by_type(NumTypes, 0);
    vector<int> str_out_by_type(NumTypes, 0);
    for(int to=0; to<names.size(); to++) {
	int str = from_to[i][to];
	if (str == 0) 
	    continue;
	if (!isdigit(names[to][0])) {
	    fprintf(logfile, "   To %15s, strength %d\n", names[to], str);
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
	    fprintf(logfile, " From %15s, strength %d\n", names[from], str);
            int ty = TypeIndex(names[from]);
            ins_by_type[ty]++;
            str_ins_by_type[ty] += str;
	    }
	}
    for(int k=0; k<NumTypes; k++) {
        if (ins_by_type[k] > 0 || out_by_type[k] > 0) {
	    fprintf(logfile, "  Type %15s, ins %d cells, total str %d,  outs %d cells, total str %d\n", CellTypes[k].name,
	     ins_by_type[k], str_ins_by_type[k], out_by_type[k], str_out_by_type[k] );
	    }
	}
    }
}
// End of completeness -----------------------------------------------------------
fprintf(logfile, "Read %d Tbars\n", tbs.size());
// Plot maps of Synapse spots
// These should look up the synapse IDS from the names, but we cheat here.
//MapXYs(tbs, 568121, "bf", name_map, names);
//MapXYs(tbs, 571372, "du", name_map, names);
//MapXYs(tbs, 597861, "fb", name_map, names);
//MapXYs(tbs, 598856, "ud", name_map, names);
// Look for some particular synapses
{
int pre = 150525;
int post=  54977;
for(int i=0; i<tbs.size(); i++) {
    if (tbs[i].body_id != pre)
        continue;
    for(int j=0; j<tbs[i].partners.size(); j++) {
        if (tbs[i].partners[j].body_id == post)
            fprintf(logfile, "TRACE %d %d %d\n", int(tbs[i].partners[j].pt.x), int(tbs[i].partners[j].pt.y), tbs[i].partners[j].z);
        }
    }
}

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
fprintf(logfile, "Strength/TBar ratio\n");
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
	    fprintf(logfile, "  From %11s, to %11s, ratio=%.3f (%4d conn, %4d Tbars)\n", names[from], names[to], double(str)/look->n, str, look->n);
	}
    }
fprintf(logfile, "End of ratio\n");
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
		fprintf(logfile, "No named connections for %s, layer M%d\n", names[T_id], med_lay);
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
fprintf(logfile, "2:Strength/TBar ratio\n");
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
	    fprintf(logfile, "What?? can't find %s %s %d %d\n", names[from], names[to], from, to);
	else {
	    fprintf(logfile, "  From %11s, to %11s, ratio=%.3f (%4d conn, %4d Tbars)\n", names[from], names[to], double(str)/look->n, str, look->n);
            string name = string(BaseName(names[from])) + ":" + string(BaseName(names[to]));
	    fprintf(logfile, "Inserting %s\n", name.c_str() );
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
fprintf(logfile, "2:End of ratio\n");
}
fprintf(logfile, "There are %d unique ID pairs\n", ss.size());
// Print a table of cell types and synapse counts
fprintf(logfile, "Cell type     M1    M2    M3    M4    M5    M6    M7    M8    M9   M10\n");
fprintf(logfile, "----------  ----  ----  ----  ----  ----  ----  ----  ----  ----  ----\n");
vector<int> HowManyTypes(11,0);
for(int i=0; i<NumTypes; i++) {
    fprintf(logfile, "%10s", CellTypes[i].name);
    for(int j=1; j<=10; j++) {
	fprintf(logfile, "%6d", s_vs_ml[i][j] );
        if (s_vs_ml[i][j] >= 10)
	    HowManyTypes[j]++;
	}
    fprintf(logfile, "\n");
    }
fprintf(logfile, "----------  ----  ----  ----  ----  ----  ----  ----  ----  ----  ----\n");
fprintf(logfile, "Types >=10");
for(int j=1; j<=10; j++)
    fprintf(logfile, "%6d", HowManyTypes[j]);
fprintf(logfile, "\n");
// Print the names of all cells with Tbars
fprintf(logfile, "\nTbar counts\n");
for(int i=0; i<names.size(); i++) {
    if (Tbar_count[i] > 0)
	fprintf(logfile, "Cell %12s has %d Tbars\n", names[i], Tbar_count[i]);
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
fprintf(logfile, "TBars for Michael\n");
for(int i=0; i<NumTypes; i++)
    fprintf(logfile, "%15s %4d\n", CellTypes[i].name, TBs[i] );
fprintf(logfile, "End Michael\n");
}
// Write a file for clustering attempts.  Compute a feature vector for each named cell.
fprintf(logfile, "Clustering test\n");
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
    fprintf(logfile, "%12s :", names[i]);
    fflush(stdout);
    for(int k=1; k<=10; k++)
	fprintf(logfile, "%5d", ins[k]);
    fprintf(logfile, "    ");
    for(int k=1; k<=10; k++)
	fprintf(logfile, "%5d", outs[k]);
    fprintf(logfile, "%8.2f %8.2f\n", volumes[i], areas[i]);
    fflush(stdout);
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
    fprintf(logfile, "Read %d overlaps\n", i);
    vector<sort_overlap> sort_over(Overlaps.size());
    int j=0;
    for(set<overlap>::iterator i = Overlaps.begin(); i != Overlaps.end(); i++)
	sort_over[j++].it = i;
    sort(sort_over.begin(), sort_over.end() );
    fprintf(logfile, "Sorted!\n");
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
	 fprintf(logfile, "%9s %9s %10.4f - %3d %3d %3d  - %5.2f %5.2f %5.2f per square micron\n", name1, name2, sort_over[i].it->area, 
	 n1, n2, sum, n1/area, n2/area, sum/area);
	}
    return 42;
    }

// Calculate an graph for simulation.
vector<bool> wanted(NumTypes, false);
vector<bool> trace (NumTypes, true);
const int threshold = 3;

int root = 568121;
CircuitBacktrace(name_map[root], from_to, names, wanted, trace, threshold);
root = 597861;
CircuitBacktrace(name_map[root], from_to, names, wanted, trace, threshold);

// sort the names by their assigned color.  (Better than names, which sort oddly in Mi1, Mi15, Mi1b, etc.)
vector<sortname> sn(M);
for(int k=0; k<names.size(); k++) {
    sn[k].index = k;
    // Normally, we sort by color.  This is an integer 0..2^24-1.  All cells of the same type have the same color,
    // so this sorts by cell type.  Next we sort by name - this puts the columns A-F at the top of the table.
    // But for non-columnar cells, columns make no sense, so we'd like to sort by size.  We can do this by
    // adding a small (less than 1.0) fraction to the color.  Since color overrides name, this will sort them
    // correctly.  And for columnar cells, color will still be identical, so sort will go to second key, name.
    sn[k].color = color(names[k]);
    fprintf(logfile, "Tweaking %s by %f\n", names[k], volumes[k]/1e6);
    sn[k].color += 1-volumes[k]/1e6;
    sn[k].name  = names[k];
    char *last = names[k] + strlen(names[k]) -1;
}
sort(sn.begin(), sn.end() );

// reverse the mapping
for(int k=0; k<names.size(); k++) {
    //              a->  b->  c->  d->  e->  f->
    char *last = names[k] + strlen(names[k]) -1;
}
vector<sort_down> sa(M);  // create a sort array
int cc = 0;            // current color

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
vector<MeanStd>  in_stats(NumTypes, MeanStd());  // stats for all connections (not just ones in table)

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
            in_stats, Tbar_count);
        if (out_row_ids.size() > 0) {
	    write_html(names, rev_map, out_row_ids, "out", out_rows, out_rows_dat, out_rows_recip, out_row_sums, out_row_vols, out_row_area, out_row_inc,
            out_stats, Tbar_count);
            // add the entries of the current table to the single column (from,to) matrix
            int ti = TypeIndex(names[out_row_ids[0]]);
	    }
	    
	in_row_ids.clear();   in_rows.clear();  in_rows_dat.clear();  in_row_sums.clear();  in_rows_recip.clear();
	out_row_ids.clear(); out_rows.clear(); out_rows_dat.clear(); out_row_sums.clear(); out_rows_recip.clear();
        in_row_area.clear(); out_row_area.clear(); in_row_vols.clear(); out_row_vols.clear();
        in_row_inc.clear(); out_row_inc.clear();
        out_stats.clear(); out_stats.resize(NumTypes, MeanStd());
        in_stats.clear();  in_stats.resize(NumTypes, MeanStd());
	fprintf(logfile, "\n");
	}
    if (ii == M)
	break;
    fprintf(logfile, "----------------- Starting analysis of '%s'\n", names[i]);
    cc = color(names[i]);
    fprintf(logfile, "----------------- Color is %d %x\n", cc,cc);

    int sos = 0;  // sum of strengths
    for(int j=0; j<M; j++) {
	sa[j].index = j;
        int str = from_to[i][j];
	sa[j].val = str;
	sos += str;
        if (str > 0) {
            int tii = TypeIndex(names[i]);
            int tij = TypeIndex(names[j]);
	    out_stats[TypeIndex(names[j])].Element(str);
            }
	}
    sort(sa.begin(), sa.end());
    if (sa[0].val > 0) {
        fprintf(logfile, "cell %11s -> ", names[i]);
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
        fprintf(logfile, "\n");
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
            }
	}
    sort(sa.begin(), sa.end());
    if (sa[0].val > 0) {
        fprintf(logfile, "cell %11s <- ", names[i]);
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
        fprintf(logfile, "\n");
	}
    }
    fclose(logfile);
    printf("Complete\n");
}
