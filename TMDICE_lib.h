#include"deps.h"
#include"globalconsts.h"


#ifndef TMDICE_lib	
#define TMDICE_lib

typedef double (*fkt)(double);
typedef double (*fkt2)(double,double);



class fourmom
{
    public:
    double px,py,pz,p0;
    fourmom(){};
    fourmom(double p00,double px0,double py0,double pz0){p0=p00;px=px0;py=py0;pz=pz0;};
    fourmom(const fourmom& org){p0=org.p0;px=org.px;py=org.py;pz=org.pz;};
    fourmom& operator=(const fourmom& org )
	{
		if(this!=&org)
        {
            p0=org.p0;px=org.px;py=org.py;pz=org.pz;
        }
        return *this;
    };
    bool operator==(const fourmom& org)
    {
        return p0==org.p0 and px==org.px and py==org.py and pz==org.pz;
    };
    fourmom operator+(const fourmom& rhs ) const
    {
        fourmom tmp ( p0+rhs.p0,px+rhs.px, py+rhs.py,pz+rhs.pz) ;
        return tmp ;
    };
    fourmom operator*(const double& a) const
    {
    fourmom tmp ( a*p0,a*px, a*py,a*pz) ;
    return tmp ;        
    };
    double operator*(const fourmom& rhs) const
    {
        return p0*rhs.p0-px*rhs.px-py*rhs.py-pz*rhs.pz;
    };
    double abs()
    {
        fourmom  tmp=*this;
        return sqrt(tmp*tmp);
    };
    double sk3(fourmom p2)
    {
        return px*p2.px+py*p2.py+pz*p2.pz;
    };
    /**@brief Pseudorapidity of jet particles. */
    double eta()
    {
        double p=sqrt(px*px+py*py+pz*pz);
        return 0.5*log((p+pz)/(p-pz));
    };

    /**@brief Rapidity of jet particles. */
    double y()
    {
        return 0.5*log((p0+pz)/(p0-pz));
    };

    /**@brief Azimutal angle of jet particles. */
    double phi()
    {
        return atan2(py,px);
    };

    /**@brief Transverse momentum of jet particles. */
    double pt()
    {
        return sqrt(px*px+py*py);
    };
};

class fourmat
{
    private:
    double a[4][4];
    public:
    fourmat(){};
    void set(int i, int j, double x){if(0<=i and i<4 and j>=0 and j<4){a[i][j]=x;} };
    fourmat& operator=(const fourmat& org )
	{
		if(this!=&org)
        {
            for(int i=0;i<=3;i++)
            {
                for(int j=0;j<=3;j++)
                {
                    a[i][j]=org.a[i][j];
                }
            }
        }
        return *this;
    };
    fourmat operator*(const fourmat&B) const
    {
        fourmat C;
        for(int i=0;i<=3;i++)
        {
            for(int j=0;j<=3;j++)
            {
                C.a[i][j]=a[i][0]*B.a[0][j]+a[i][1]*B.a[1][j]+a[i][2]*B.a[2][j]+a[i][3]*B.a[3][j];
            }
        }
        return C;
    };
    fourmom operator*(const fourmom&v)const
    {
        fourmom w;
        w.p0=a[0][0]*v.p0+a[0][1]*v.px+a[0][2]*v.py+a[0][3]*v.pz;
        w.px=a[1][0]*v.p0+a[1][1]*v.px+a[1][2]*v.py+a[1][3]*v.pz;
        w.py=a[2][0]*v.p0+a[2][1]*v.px+a[2][2]*v.py+a[2][3]*v.pz;
        w.pz=a[3][0]*v.p0+a[3][1]*v.px+a[3][2]*v.py+a[3][3]*v.pz;
        return w;
    }
    fourmat transpose()
    {
        fourmat B;
        for(int i=0; i<=3; i++)
        {
            for(int j=0;j<=3;j++)
            {
                B.set(i,j,a[j][i]);
            }
        }
        return B;
    };
    double get(int i, int j){return a[i][j];};
};

const double C_F=4.0f/3.0f;
const double T_R=0.5f;
const double C_A=3.0f;


extern int n;
extern double gaussx[10000],gaussw[10000];
extern double N_F,N_c;

void fillgauss(double (&gaussx)[10000],double (&gaussw)[10000]);

template<typename T, typename ...Rest> double integral(  double f(  T x,Rest...rest), double a,  double b,T x, Rest...rest );

double rand01(void);


template<typename T, typename ...Rest> void fill_partitionfunction_v2(  double f(  T x,Rest...rest), double vmin,  
double vmax, double vanz,std::map<double,double> &Finv, double &Fmax,T x, Rest...rest );

std::unordered_map<int,double> discretize_grid(std::map<double,double>oldmap,double mapmax, int anz);

void fill_partitionfunction_v2( fkt f, double vmin,double vmax,double vanz,std::map<double,double> &Finv, std::map<double,double>&F,double &Fmax);

void fill_partitionfunction_v2( fkt f, double vmin,double vmax,double vanz,std::map<double,double> &Finv, double &Fmax);



double select(fkt Finv, double Fmax,double wmin,double wmax);

double T(double x, int n);

std::vector<double> det_boundaries(fkt f, double y, double zmin,double zmax);

void make_c(std::map<double,double> f,double (&c)[M],double a, double b,double eee);

void make_3c(std::map<double,double>f,double (&c1)[M],double (&c2)[M],double (&c3)[M],double a,double b1,double b2,double b);

double approx_f(double x, double a, double b,double (&c)[M],double eee);

double approx_f(double x,double *c1,double *c2,double *c3,double a,double b1,double b2,double b,int M);

double interpolate_integer_inversegrid_1d(std::unordered_map<int,double> *phinvtmp,double xminus, double xplus, double r);

double select(fkt2 Finv, double Fmax,double wmin,double wmax,double t);

std::vector<double> det_boundaries(fkt f, double y, double zmin,double zmax);

std::vector<double> det_boundaries(fkt2 f, double y, double zmin,double zmax,double t);

void make_3c(std::map<double,double>f,double (&c1)[M],double (&c2)[M],double (&c3)[M],double a,double b1,double b2,double b);


void cheb_coeff_split(std::map<double,std::map<double,fkt2>>kernel,std::map<double,std::map<double,std::map<double,double>>>ph_max,
std::map<double,std::map<double,std::map<double,std::map<double,double>>>>phi_inv,int f1, int f2
, double ymin,double (&c1)[3][3][tanzmax][M],double (&c2)[3][3][tanzmax][M],double (&c3)[3][3][tanzmax][M]
,std::map<double,std::map<double,std::map<double,double>>> (&b1),std::map<double,std::map<double,std::map<double,double>>> (&b2),
std::map<double,std::map<double,std::map<double,double>>> (&b),double t);

double interpolate_cdf(double z, double z_max, std::unordered_map<int,double> cdf, int cdfanz,double cdfmax=1. ,double z_min=0.);

void print_cdfs(double f1,double f2);

void GetTransverseDirections(fourmom p1,fourmom &p2, fourmom &p3);

fourmom Rotate(fourmom axis,fourmom vec);//vector is the vector that will be rotated;

#include"TMDICE_lib.tpp"
#endif

