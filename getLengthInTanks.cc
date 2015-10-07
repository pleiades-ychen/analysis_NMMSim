// getLengthInTanks.cc
// Author: Yu Chen (Feb. 2013)
// Plot distribution of lengths that muons would travel inside water tanks

// //	to compile, do:
//	g++ -O2 `$ROOTSYS/bin/root-config --libs` -I$ROOTSYS/include getLengthInTanks.cc -o getLengthInTanks
//	to run, there is one arugments:
//	./getLengthInTanks <MuonDataListFile>

#define _CRTDBG_MAP_ALLOC
#include <stdlib.h>
//#include <crtdbg.h>
//  ROOT includes
//
#include "TChain.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TLeaf.h"
#include "TText.h"
#include "TH1.h"
#include "TH2.h"
#include "TRandom3.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TLegend.h"
#include "TVector3.h"
//
//  C/C++ includes
//
#include <fstream>
#include <stdio.h>
#include <math.h>
#include <string>
#include <iostream>
#include <sstream>
//
//  Definitions
//
#define DEBUGGING 0

using namespace std;

class TrackLine {
    double x0,y0,z0; // primary position
    double dx,dy,dz; // direction vector
    
public:
    TrackLine(double,double,double,double,double,double);
    TVector3 GetPoint(double,char);
};

TrackLine::TrackLine(double posX, double posY, double posZ, double dirX, double dirY, double dirZ){
    x0 = posX;
    y0 = posY;
    z0 = posZ;
    dx = dirX;
    dy = dirY;
    dz = dirZ;
}

TVector3 TrackLine::GetPoint(double coord,char coordName) {
    double x,y,z;
    double ratio;
    if (coordName=='x') {
        x = coord;
        ratio = (x - x0)/dx;
        y = ratio*dy + y0;
        z = ratio*dz + z0;
    } else if(coordName=='y') {
        y = coord;
        ratio = (y - y0)/dy;
        x = ratio*dx + x0;
        z = ratio*dz + z0;
    } else if(coordName=='z'){
        z = coord;
        ratio = (z - z0)/dz;
        x = ratio*dx + x0;
        y = ratio*dy + y0;
    } else {
        cout<<"Please provide coordinate name. Candidates: 'x', 'y', or 'z'.\n";
        return TVector3(x0,y0,z0);
    }
    
    return TVector3(x,y,z);
}
/////////////////////////////////////////////////////////////////////
class Plane {
    int norm;
    double para;
    
public:
    Plane();
    Plane(double,int);
    int Norm(){return (norm);}
    char NormName();
    double Para(){return (para);}
};

Plane::Plane(){ }

Plane::Plane(double coord,int direction){
    norm = direction;
    para = coord;
}

char Plane::NormName(){
    if (norm==0) return 'x';
    if (norm==1) return 'y';
    if (norm==2) return 'z';
}
/////////////////////////////////////////////////////////////////////
class Rectangle {
    double center[2];
    double side[2];
    
public:
    Rectangle();
    Rectangle(double,double,double,double);
    // i = 0,1; for the two dimentions of Rectangle
    // Bottom(i) and Top(i) mean Lower and Upper limits of the rectangle's ranges for each dimensions.
    bool Inside(double,double);
    double Bottom(int i){return (center[i]-0.5*side[i]);}
    double Top(int i){return (center[i]+0.5*side[i]);}
};

Rectangle::Rectangle(){ }

Rectangle::Rectangle(double cent1,double cent2, double w, double h){
    center[0] = cent1;
    center[1] = cent2;
    side[0] = w;
    side[1] = h;

}

bool Rectangle::Inside(double a, double b){
    if ( (Bottom(0)<=a && a<=Top(0)) && (Bottom(1)<=b && b<=Top(1)) ) {
        return true;
    } else{
        return false;
    }
}
///////////////////////////////////////////////////////////////////////
class RightBox {
    double center[3];
    double side[3];
public:
    RightBox(double,double,double,double,double,double);
    Rectangle GetFace(int);
    Plane GetPlane(int);
};

RightBox::RightBox(double cent1, double cent2, double cent3, double a, double b, double c) {
    center[0] = cent1;
    center[1] = cent2;
    center[2] = cent3;
    side[0] = a;
    side[1] = b;
    side[2] =c;
}

Rectangle RightBox::GetFace(int i){
    // Only need 3 faces, identify the front and back the same, the left and right the same, and the top and bottom the same. 
    // i=0,1,2 for x,y,z
    int dimA = (i+1) % 3;
    int dimB = (i+2) % 3;
    return Rectangle(center[dimA],center[dimB],side[dimA],side[dimB]);
}

Plane RightBox::GetPlane(int i){
    // There're 6 planes associated with a box
    // i =0,1,...,5 for front, back (x), right, left (y), top, and bottom (z)
    int dim = i/2; // dim = 0,1,2 for x,y,z
    int sign = 1 - (i % 2)*2;
    double coord = center[dim] + sign * side[dim]*0.5;
    return Plane(coord,dim);
}
///////////////////////////////////////////////////////////////////////////

int main(int argc, char** argv){
    
    int event , type;
    double energy, posX, posY, posZ, dirX, dirY, dirZ;
    
    TH1F* h = new TH1F("hPassLength"
                       ,"Distribution of Pass Lengths of Muon Tracks in Tanks;Pass Length (cm);Counts/cm"
                       ,300,0,300); //cm
    const double offsetX = 594.;
    const double offsetY = 198.4;
    const double offsetZ = -287.72;
    RightBox * northTank = new RightBox(offsetX+ 0.,offsetY+ 63.5, offsetZ+ 1.27
                                        ,243.84,121.92,78.105); // cm
    RightBox * southTank = new RightBox(offsetX+ 0., offsetY+ -63.5, offsetZ+ 6.35
                                        ,243.84,121.92,78.105);
    // take care maybe x, y order is wrong
    Plane ntp[6],stp[6];
    for (int i=0;i<6;i++){
        ntp[i] = northTank->GetPlane(i);
        stp[i] = southTank->GetPlane(i);
    }
    Rectangle ntf[3],stf[3];
    for (int i=0; i<3; i++) {
        ntf[i] = northTank->GetFace(i);
        stf[i] = southTank->GetFace(i);
    }
    
    int lastevent = -1;
    ifstream InList;
    InList.open(argv[1]);
    while (!InList.eof()){
        char filename[100];
        InList>>filename;
        ifstream muonfile;
        muonfile.open(filename);
        while (!muonfile.eof()){
            muonfile >> event >> type >> energy >> posX >> posY >> posZ >> dirX >> dirY >> dirZ;
            // take care
            if(event==lastevent) continue;
            
            double PassLength,passNT,passST;
            TrackLine * muon = new TrackLine(posX,posY,posZ,dirX,dirY,dirZ);
            
            vector<TVector3> xingNT,xingST;
            // north tank
            for (int i=0; i<6; i++) {
                TVector3 xingPl = muon->GetPoint( ntp[i].Para() , ntp[i].NormName());
                if ( ntf[i/2].Inside(xingPl( (i/2+1)%3 ) , xingPl( (i/2+2)%3 ) ) ==true) {
                    int NumSame=0;
                    for (int j=0; j<xingNT.size(); j++) {
                        if (xingPl==xingNT[j]) NumSame++;
                    }
                    if (NumSame==0) xingNT.push_back(xingPl);
                }
            }
            
            if (xingNT.size()==2) {
                TVector3 passVecN = xingNT[0]-xingNT[1];
                passNT = passVecN.Mag();
            } else if (xingNT.size()==0){
                passNT = 0.;
            } else {
                cout<<"Error on crossing judgment on north tank. You get to debug...\n";
            }
            
            // south tank
            for (int i=0; i<6; i++) {
                TVector3 xingPl = muon->GetPoint( stp[i].Para() , stp[i].NormName());
                if ( stf[i/2].Inside(xingPl( (i/2+1)%3 ) , xingPl( (i/2+2)%3 ) ) ==true) {
                    int NumSame=0;
                    for (int j=0; j<xingST.size(); j++) {
                        if (xingPl==xingST[j]) NumSame++;
                    }
                    if (NumSame==0) xingST.push_back(xingPl);
                }
            }
            if (xingST.size()==2) {
                TVector3 passVecS = xingST[0]-xingST[1];
                passST = passVecS.Mag();
            } else if (xingST.size()==0){
                passST = 0.;
            } else {
                cout<<"Error on crossing judgment on south tank. You get to debug...\n";
            }
            
            PassLength = passNT + passST;
            
            h->Fill(PassLength);
            
            lastevent = event;
        }// An Event
        
    }// A File
    
    TFile * OutFile = (TFile*) gROOT->FindObject("histo.root");
    if(OutFile) OutFile->Close();
	OutFile = new TFile("histo.root","RECREATE");
    
    TCanvas* c1 = new TCanvas("c1");
    h->Draw();
    h->Write();
    OutFile->Close();
}





























