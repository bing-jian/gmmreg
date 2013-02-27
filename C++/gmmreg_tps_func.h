/*=========================================================================
$Author: bing.jian@gmail.com $
$Date: 2010-05-10 01:32:09 -0400 (Mon, 10 May 2010) $
$Revision: 127 $
=========================================================================*/

/**
 * \file gmmreg_tps_func.h
 * \brief  The declaration of the class gmmreg_tps_func
 */


#ifndef gmmreg_tps_func_h
#define gmmreg_tps_func_h


#include <vnl/vnl_cost_function.h>
#include "gmmreg_base.h"

/**
 * \class gmmreg_tps_func
 * \brief  the abstract base class for the L2/KC-based cost function
 */


class gmmreg_tps_func : public vnl_cost_function{

public:
    gmmreg_tps_func(): vnl_cost_function() {}

    virtual double eval(double&f1, double &f2, vnl_matrix<double> &g1, vnl_matrix<double> &g2) = 0;
    double f(const vnl_vector<double>& x);
    void gradf(const vnl_vector<double>& x, vnl_vector<double>& g);

    gmmreg_base* gmmreg;
    inline void set_gmmreg(gmmreg_base* gmmreg){
        this->gmmreg = gmmreg;
        this->m = gmmreg->m;
        this->n = gmmreg->n;
        this->d = gmmreg->d;
        gradient1.set_size(m,d);
        gradient2.set_size(m,d);
    }
    inline void set_scale(double scale){this->scale = scale;}
    inline double get_scale(){return this->scale;}
    inline void set_lambda(double lambda){this->lambda = lambda;}
    inline double get_lambda(){return this->lambda;}

    void prepare_param_gradient();
    bool fix_affine;
    inline void set_fix_affine(bool fix_affine){
        this->fix_affine = fix_affine;
        if (fix_affine) {
            set_number_of_unknowns((n-d-1)*d); //dim = (n-d-1)*d;
        }
        else  {
            set_number_of_unknowns(n*d); //dim = n*d;
        }
        //gmmreg->dim = get_number_of_unknowns(); //dim;
    }
    inline bool get_fix_affine(){return this->fix_affine;}

    virtual ~gmmreg_tps_func() {}

protected:
    vnl_matrix<double> gradient;

private:
    double scale, lambda;
    int m,n,d;
    vnl_matrix<double> gradient1, gradient2, grad_all;
    

};

class gmmreg_tps_L2_func : public gmmreg_tps_func{
    double eval(double&f1, double &f2, vnl_matrix<double> &g1, vnl_matrix<double> &g2);
};


class gmmreg_tps_KC_func : public gmmreg_tps_func{
    double eval(double&f1, double &f2, vnl_matrix<double> &g1, vnl_matrix<double> &g2);
};

#endif //#ifndef gmmreg_tps_func_h_

