// calcLonLat2.cpp : Defines the entry point for the console application.
//这个程序由已知惯性系的原点和目标在惯性系下的坐标计算目标的ＧＰＳ经纬度

#include "stdafx.h"
#include "stdio.h"
#include "math.h"

double a=6378137,b=6356752.3142;
double af=a*a,bf=b*b;
//double e2f=af/bf-1;
double e1f=1-bf/af;	//0.00669438第一篇心率的平方

double long_missle;
double lat_missle;
double long0;
double lat0;
double RRn[3];
double Cen[3][3];
double R0e[3];
double RRe[3];
double RR_fashe[3];
double ECEFX,ECEFY,ECEFZ;
double Lat_ld,Lon_ld,Alt_ld;
double X_ld,Y_ld,Z_ld;
double PosX0,PosY0,PosZ0;
double RX_ld,RY_ld,RZ_ld;
double X_m,Y_m,Z_m;
double RX_m,RY_m,RZ_m;
double dzE,dxN;
#define PI 3.14159265358979
#define deg_to_rad (PI/180.0)
#define rad_to_deg (180.0/PI)

#define RADA 6378137.0
#define RADB 6356752.3142

/*/大地坐标系到ECEF坐标的转换//GPSTo84
void llhtoECEF(float lat,float lon,float h)
{
	float af,bf,ef,clat,slat,NN;

	af = RADA*RADA;
	bf = RADB*RADB;
	ef = 1-bf/af;
	clat = cos(lat);
	slat = sin(lat);
	NN = RADA/sqrt(1-ef*slat*slat);

	//输出ECEF坐标
	ECEFX = (NN+h)*clat*cos(lon);
	ECEFY = (NN+h)*clat*sin(lon);
	ECEFZ = (NN*(1.0-ef)+h)*slat;
}
*/


int main(int argc, char* argv[])
{
	int i;
	double N,h;
	double slat0,slong0;
	double clat0,clong0;
	double slat,clat;
	double lat1,dat2,lat3;
	double long1,long2,long3;
	ECEFX = 0.0;			//ECEF--地球固联坐标系
	ECEFY = 0.0;
	ECEFZ = 0.0;

	//目标位置的惯性系坐标
	RRn[0] = 529.77;
	RRn[1] = 0;
	RRn[2] = 281.683;
	//原点的经纬度
	lat0 = 40.639722;
	long0 = 100.4805556;
	//真实落点坐标
	//Lat_ld = 40;
	//Lon_ld = 100;
	//Alt_ld= 22.2;

	//计算惯性系和大地坐标系的转换矩阵
	//这个矩阵是从地心到地理的，如果要从地理到地心需要取矩阵的逆
	//但是对于转换矩阵这个正交矩阵来说，逆就是它的转置，因此可以
	//看见计算式中算式实际是取得这个矩阵的逆
	slat0 = sin(lat0/rad_to_deg);
	clat0 = cos(lat0/rad_to_deg);
	slong0 = sin(long0/rad_to_deg);
	clong0 = cos(long0/rad_to_deg);

	Cen[0][0] = -slat0*clong0;
	Cen[0][1] = -slat0*slong0;
	Cen[0][2] = clat0;
	Cen[1][0] = clat0*clong0;
	Cen[1][1] = clat0*slong0;
	Cen[1][2] = slat0;
	Cen[2][0] = -slong0;
	Cen[2][1] = clong0;
	Cen[2][2] = 0;

	//计算已知点（原点）在大地坐标系中的坐标R0e
	N=a/sqrt(1-e1f*slat0*slat0);
	R0e[0] = N*clat0*clong0;
	R0e[1] = N*clat0*slong0;
	R0e[2] = N*(1-e1f)*slat0;
	//计算要计算的点在大地坐标系的分量RRe//该步骤实际是将惯性系（地理系）点转换为地心系后加上原点坐标
	RRe[0] = R0e[0]+(Cen[0][0]*RRn[0]+Cen[1][0]*RRn[1]+Cen[2][0]*RRn[2]);
	RRe[1] = R0e[1]+(Cen[0][1]*RRn[0]+Cen[1][1]*RRn[1]+Cen[2][1]*RRn[2]);
	RRe[2] = R0e[2]+(Cen[0][2]*RRn[0]+Cen[1][2]*RRn[1]+Cen[2][2]*RRn[2]);

	//根据RRe计算经纬度 84ToＧＰＳ
	double x = RRe[0],y=RRe[1],z=RRe[2];
	long_missle=atan2(y,x);		//计算经度

	//循环计算纬度
	double p = sqrt(x*x+y*y);
	N=a,h=12800;
	for(i=0;i<8;i++)
	{
		lat_missle = atan(z/p/(1-e1f*(N/(N+h))));
		clat = cos(lat_missle);
		h = p/clat - N;
		slat = sin(lat_missle);
		N = a/sqrt(1-e1f*slat*slat);
	}

	//计算结果输出
	long_missle = long_missle*rad_to_deg;
	lat_missle = lat_missle*rad_to_deg;
	//把经纬度转换成度分秒
	int long_missle1,long_missle2;
	int lat_missle1,lat_missle2;
	double long_missle3,lat_missle3,temp1;

	long_missle1 = floor(long_missle);
	temp1 = (long_missle-long_missle1)*60;
	long_missle2=floor(temp1);
	long_missle3 = (temp1-long_missle2)*60;

	lat_missle1 = floor(lat_missle);
	temp1=(lat_missle-lat_missle1)*60;
	lat_missle2 = floor(temp1);
	lat_missle3 = (temp1-lat_missle2)*60;

//	llhtoECEF(Lat_ld/rad_to_deg,Lon_ld/rad_to_deg,Alt_ld);
//	X_ld = ECEFX;
//	Y_ld = ECEFY;
//	Z_ld = ECEFZ;
//	PosX0=R0e[0];
//	PosY0=R0e[1];
//	PosZ0=R0e[2];

//	RX_ld=Cen[0][0]*(X_ld-PosX0)+Cen[0][1]*(Y_ld-PosY0)+Cen[0][2]*(Z_ld-PosZ0);
//	RY_ld=Cen[1][0]*(X_ld-PosX0)+Cen[1][1]*(Y_ld-PosY0)+Cen[1][2]*(Z_ld-PosZ0);
//	RZ_ld=Cen[2][0]*(X_ld-PosX0)+Cen[2][1]*(Y_ld-PosY0)+Cen[2][2]*(Z_ld-PosZ0);
	//遥测落点坐标
//	llhtoECEF(lat_missle/rad_to_deg,long_missle/rad_to_deg,Alt_ld);
//	X_m = ECEFX;
//	Y_m = ECEFY;
//	Z_m = ECEFZ;
//	RX_m=Cen[0][0]*(X_m-PosX0)+Cen[0][1]*(Y_m-PosY0)+Cen[0][2]*(Z_m-PosZ0);
//	RY_m=Cen[1][0]*(X_m-PosX0)+Cen[1][1]*(Y_m-PosY0)+Cen[1][2]*(Z_m-PosZ0);
//	RZ_m=Cen[2][0]*(X_m-PosX0)+Cen[2][1]*(Y_m-PosY0)+Cen[2][2]*(Z_m-PosZ0);

	//惯性系两点间差距
//	dzE=RRn[2]-RZ_ld;
//	dxN=RRn[0]-RX_ld;

	printf("落点纬度=%lf度\n",lat_missle);
	printf("落点经度=%lf度\n\n\n",long_missle);

	printf("落点纬度=%d度  %d分  %d秒\n",lat_missle1,lat_missle2,(int)lat_missle3);
	printf("落点经度=%d度  %d分  %d秒\n\n\n",long_missle1,long_missle2,(int)long_missle3);

//	printf("北向相差距离=%d",dxN);
//	printf("东向相差距离=%d",dzE);
	return 0;
}
