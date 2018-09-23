// calcLonLat2.cpp : Defines the entry point for the console application.
//�����������֪����ϵ��ԭ���Ŀ���ڹ���ϵ�µ��������Ŀ��ģǣУӾ�γ��

#include "stdafx.h"
#include "stdio.h"
#include "math.h"

double a=6378137,b=6356752.3142;
double af=a*a,bf=b*b;
//double e2f=af/bf-1;
double e1f=1-bf/af;	//0.00669438��һƪ���ʵ�ƽ��

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

/*/�������ϵ��ECEF�����ת��//GPSTo84
void llhtoECEF(float lat,float lon,float h)
{
	float af,bf,ef,clat,slat,NN;

	af = RADA*RADA;
	bf = RADB*RADB;
	ef = 1-bf/af;
	clat = cos(lat);
	slat = sin(lat);
	NN = RADA/sqrt(1-ef*slat*slat);

	//���ECEF����
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
	ECEFX = 0.0;			//ECEF--�����������ϵ
	ECEFY = 0.0;
	ECEFZ = 0.0;

	//Ŀ��λ�õĹ���ϵ����
	RRn[0] = 529.77;
	RRn[1] = 0;
	RRn[2] = 281.683;
	//ԭ��ľ�γ��
	lat0 = 40.639722;
	long0 = 100.4805556;
	//��ʵ�������
	//Lat_ld = 40;
	//Lon_ld = 100;
	//Alt_ld= 22.2;

	//�������ϵ�ʹ������ϵ��ת������
	//��������Ǵӵ��ĵ�����ģ����Ҫ�ӵ���������Ҫȡ�������
	//���Ƕ���ת�������������������˵�����������ת�ã���˿���
	//��������ʽ����ʽʵ����ȡ������������
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

	//������֪�㣨ԭ�㣩�ڴ������ϵ�е�����R0e
	N=a/sqrt(1-e1f*slat0*slat0);
	R0e[0] = N*clat0*clong0;
	R0e[1] = N*clat0*slong0;
	R0e[2] = N*(1-e1f)*slat0;
	//����Ҫ����ĵ��ڴ������ϵ�ķ���RRe//�ò���ʵ���ǽ�����ϵ������ϵ����ת��Ϊ����ϵ�����ԭ������
	RRe[0] = R0e[0]+(Cen[0][0]*RRn[0]+Cen[1][0]*RRn[1]+Cen[2][0]*RRn[2]);
	RRe[1] = R0e[1]+(Cen[0][1]*RRn[0]+Cen[1][1]*RRn[1]+Cen[2][1]*RRn[2]);
	RRe[2] = R0e[2]+(Cen[0][2]*RRn[0]+Cen[1][2]*RRn[1]+Cen[2][2]*RRn[2]);

	//����RRe���㾭γ�� 84To�ǣУ�
	double x = RRe[0],y=RRe[1],z=RRe[2];
	long_missle=atan2(y,x);		//���㾭��

	//ѭ������γ��
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

	//���������
	long_missle = long_missle*rad_to_deg;
	lat_missle = lat_missle*rad_to_deg;
	//�Ѿ�γ��ת���ɶȷ���
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
	//ң���������
//	llhtoECEF(lat_missle/rad_to_deg,long_missle/rad_to_deg,Alt_ld);
//	X_m = ECEFX;
//	Y_m = ECEFY;
//	Z_m = ECEFZ;
//	RX_m=Cen[0][0]*(X_m-PosX0)+Cen[0][1]*(Y_m-PosY0)+Cen[0][2]*(Z_m-PosZ0);
//	RY_m=Cen[1][0]*(X_m-PosX0)+Cen[1][1]*(Y_m-PosY0)+Cen[1][2]*(Z_m-PosZ0);
//	RZ_m=Cen[2][0]*(X_m-PosX0)+Cen[2][1]*(Y_m-PosY0)+Cen[2][2]*(Z_m-PosZ0);

	//����ϵ�������
//	dzE=RRn[2]-RZ_ld;
//	dxN=RRn[0]-RX_ld;

	printf("�������γ��=%lf��\n",lat_missle);
	printf("������㾭��=%lf��\n\n\n",long_missle);

	printf("�������γ��=%d��  %d��  %d��\n",lat_missle1,lat_missle2,(int)lat_missle3);
	printf("������㾭��=%d��  %d��  %d��\n\n\n",long_missle1,long_missle2,(int)long_missle3);

//	printf("����������=%d",dxN);
//	printf("����������=%d",dzE);
	return 0;
}
