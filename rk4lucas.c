#include <stdio.h>
#include <stdlib.h>

typedef struct tripla
{
    int tam;
    double *v;
    double *ex;
} tripla;


tripla** runge_kutta4(double x0,double y0,double t0,double t_final,double dt_max, double mi);

double f(double t);

double f1(double x, double y, double mi, double t);

double f2(double x, double y, double mi, double t);

double exata(double t);

int main(){
    int i;
    double t0,x0,y0,dt_max,t_final,t;
    double *y;
    double mi=0.07, epsilon = -1e-1;
    tripla **pp_tripla;
    t0 = 0.0;
    t_final = 100.0;
    dt_max = 0.1;
    x0 = 0.35;//(mi*mi+1)/(mi+2)+epsilon;
    y0 = 0.29;//(mi*mi+1)*(1-2*mi)/((mi+2)*(mi+2))+epsilon;
    pp_tripla = runge_kutta4(x0,y0,t0,t_final,dt_max,mi);
    printf("[");
    for(i=0;i<pp_tripla[0]->tam;i++){
        printf(" (%lf, %lf)",pp_tripla[0]->v[i],pp_tripla[1]->v[i]);
    }
    printf("]\n");
    free(pp_tripla[0]->v);
    free(pp_tripla[0]->ex);
    free(pp_tripla[1]->v);
    free(pp_tripla[1]->ex);
    free(pp_tripla[0]);
    free(pp_tripla[1]);
    free(pp_tripla);
    return 0;
}

tripla** runge_kutta4(double x0,double y0,double t0,double t_final,double dt_max,double mi){
    double k1,k2,s1,s2,l1,l2,p1,p2;
    double xaux,yaux;
    double t,dt;
    tripla ** pp_tripla;
    tripla *p_x,*p_y;
    pp_tripla = malloc(2*sizeof(tripla*)); //pp_tripla[0] = p_x e pp_tripla[1] = p_y
    if(pp_tripla == NULL)
        exit(1);
    int i,n;
    n = (int)((t_final-t0) / dt_max);
    if(t0+ n*dt_max != t_final){
        dt = (t_final-t0)/(n+1);
        n += 2;
    }
    else{
        dt = dt_max;
        n++;
    }


    p_x = pp_tripla[0] = malloc(sizeof(tripla));
    p_y = pp_tripla[1] = malloc(sizeof(tripla));
    if(p_x == NULL || p_y == NULL)
        exit(1);

    p_x->tam=n; // quantidade de pontos totais (inicial + internos + final)
    p_y->tam=n;
    p_x->v=malloc(n*sizeof(double));
    p_x->ex = malloc(n*sizeof(double));
    if( p_x->v == NULL || p_x->ex == NULL)
        exit(1);
        
    p_y->v=malloc(n*sizeof(double));
    p_y->ex = malloc(n*sizeof(double));
    if( p_y->v == NULL || p_y->ex == NULL)
        exit(1);

    t=t0;
    p_x->v[0] = x0;
    p_x->ex[0] = exata(t0);

    p_y->v[0] = y0;
    p_y->ex[0] = exata(t0);
    for(i=0; i<p_x->tam-1;i++){
        k1=dt*f1(p_x->v[i],p_y->v[i],mi,t);
        k2=dt*f2(p_x->v[i],p_y->v[i],mi,t);

        xaux = p_x->v[i]+0.5*k1*dt;
        yaux = p_y->v[i]+0.5*k2*dt;

        s1=dt*f1(xaux,yaux,mi,t+0.5*dt);
        s2=dt*f2(xaux,yaux,mi,t+0.5*dt);
        
        xaux = p_x->v[i]+0.5*s1*dt;
        yaux = p_y->v[i]+0.5*s2*dt;
        
        l1=dt*f1(xaux,yaux,mi,t+0.5*dt);
        l2=dt*f2(xaux,yaux,mi,t+0.5*dt);

        xaux = p_x->v[i]+l1;
        yaux = p_y->v[i]+l2;

        p1=dt*f1(xaux,yaux,mi,t+dt);
        p2=dt*f2(xaux,yaux,mi,t+dt);
        p_x->v[i+1] = p_x->v[i]+(k1+2*s1+2*l1+p1)*dt/6.0;
        p_y->v[i+1] = p_y->v[i]+(k2+2*s2+2*l2+p2)*dt/6.0;
        t+=dt;
        p_x->ex[i+1] = exata(t);
        p_y->ex[i+1] = exata(t);
    }

    return pp_tripla;
}

double f(double t){
    return t;
}

double f1(double x, double y, double mi, double t){
    return mi*x + y -x*x;
}

double f2(double x, double y, double mi, double t){
    return -x + mi*y +2.0*x*x;
}

double exata(double t){
    return 0.5*t*t+1.0;
}

