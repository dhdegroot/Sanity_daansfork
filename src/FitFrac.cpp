#include <FitFrac.h>

double fitfrac(double *f, double *n_c, double n, double &v, int C, double *N_c, double a, double b){
    double q,qmin,qmax,nor,normin,normax;
    int i;
    int counter(0);
	double *Q = new double [C];
	double beta = 1.0/((n+a)*v);
	double logbeta = log(beta);
	q=0;
    cerr << "I'm here, initial q " << q << endl;
	for(i=0;i<C;++i){
        Q[i] = log( N_c[i]) + n_c[i]*v -logbeta;  // TODO: Remove logbeta?
		q += N_c[i];  // TODO: This must be gene-independent, right?
    }
    cerr << "I'm here, next q" << q << endl;

	q =log(q+b);//initial guess for q. // TODO: Since b is fixed, this initial guess is always the same
    cerr << "I'm here, b " << b << ", beta" << beta << endl;
    cerr << "I'm here, initial guess " << q << endl;
    double dq = log(2.0);

    /**** get initial nor ****/
    cerr << "I'm here 10a3 " << endl;
    cerr << "Q, C, q, b: " << Q[0] << ", " << C << ", " << q << ", " << b << ", " << endl;
    nor = beta*normalization(Q,C,q)+b*exp(-q);
    cerr << "nor, beta*normalization: " << nor << ", " << beta*normalization(Q,C,q) << endl;

    cerr << "I'm here 10c1 " << nor << endl;
    exit(0);
    /**** initial c too small. Set as lower bound and increase by factor 2 until a value too large is found****/
    normin = nor;
    normax = nor;
    if(nor < 1){
        while(nor < 1){
            qmax = q;
            q = q-dq;
            nor = beta*normalization(Q,C,q)+b*exp(-q);
            cerr << "I'm here q " << q << ", normin " << nor << endl;
        }
        qmin = q;
        normin = nor;
    }
    /**** initial c too large. Set as upper bound and decrease by factor 2 until a value too small is found****/
    else{
        while(nor > 1){
            qmin = q;
            q = q+dq;
            nor = beta*normalization(Q,C,q)+b*exp(-q);
            cerr << "I'm here q " << q << ", normax " << nor << endl;
        }
        qmax = q;
        normax = nor;
    }

    cerr << "I'm here 10a31 " << endl;
        cerr << "I'm here test " << qmin << endl;
    /***now do the bisection*****/
    double tol = 0.0000001;
    double diff = 1.0;
    cerr << "I'm here 10b1 " << qmin << endl;
    cerr << "I'm here 10b2 " << qmax << endl;
    while(diff > tol){
        counter += 1;
        cerr << "I'm here 10b3 " << diff << endl;
        q = (qmax+qmin)/2;  // TODO: This update rule can be much more efficient
        cerr << "qmin, qmax, normin, normax: " << qmin << "," << qmax << "," << normin << "," << normax << endl;
    	// q = qmin + ((normin - 1) / (normin - normax)) * (qmax - qmin);
        cerr << "q:" << q << endl;
        nor = beta*normalization(Q,C,q)+b*exp(-q);
        cerr << "I'm here 10b4 " << nor << endl;
        if(nor > 1){
            qmin = q;
            normin = nor;
        }else{
            qmax = q;
            normax = nor;
        }
        diff = fabs(nor-1.0);
        if (counter >=3){exit(0);}
    }
    cerr << "I'm here 10a4 " << endl;

    cerr << "Number of iterations done: " << counter << endl;
    // q = (qmax+qmin)/2;
    q = qmin + (normin / (normin - normax)) * (qmax - qmin);
	for(i=0;i<C;++i){
        double x = Q[i] -q;
        if(x > 50){
            f[i] = beta*(x-log(x)+log(x)/x);
        }else{
            //f[i] = beta*gsl_sf_lambert_W0(exp(x));
			f[i] = beta*Fukushima::LambertW(0,exp(x));
        }
    }
	delete[] Q;
    exit(0);
    return q;
}

double normalization(double *Q,int C,double q){
    double nor = 0;
    double x;
    int i;
    for(i=0;i<C;++i){
        x = Q[i]-q;
		if(isnan(x)){
			cout << "x = nan !\n";
		}
        if(x > 50){
            nor += (x-log(x)+log(x)/x);
        }else{
            //nor += gsl_sf_lambert_W0(exp(x));
			nor += Fukushima::LambertW(0,exp(x));
        }
        if (nor > 10000){
            cerr << "i, nor: " << i << ", " << nor << endl;
            exit(0);
        }
    }
    return nor;
}

