#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>

/* function declaration */
// int secular(float d, float v, const float p, int mode);

double * secular(float *d, float *v,float p, int mode, int l_d)
{
    double * sol = (double *)malloc(sizeof(double) * l_d);
    // const float p = 1.0e-10; // define alpha value
    // int mode = 4;
    // float d[1][4] = {{4, 50, 600, 7000}}; // declare the d values
    // float v[4] = {1, 2, 3, 9};            // declare the v_vect values
    float lamda, dk, dk1, f_, f_y, g_, g_x, Tr, c1, a1, b1, init_y, Tr_k, Tr_k_1, f__, f_y_, o_, o_y_, a2, b2, c2, lamda_new, t, u_, u_y_;
    float Tao, n_cor;
    
    int k, i;
    

    // clock_t begin = clock();

    // struct timeval  tv1, tv2;
    // gettimeofday(&tv1, NULL);

    for (int i = 1; i < mode + 1; i = i + 1)
    {

        k = i;
        dk = d[k - 1];

        if (k < mode)
        {
            dk1 = d[k];
            lamda = (dk + dk1) / 2;

            /* To check if f_y is positive or negative*/
            f_ = 0;
            for (int q = 0; q < mode; q = q + 1)
            {
                f_ = f_ + v[q] / (d[q] - lamda);
            }
            f_y = f_ + p;

            /* Now, applying equation 40 and 41, choosing our initial guess */

            g_x = f_y - v[k - 1] / (d[k - 1] - lamda);
            Tr = d[k] - d[k - 1];
            c1 = g_x;

            /* first a,b and c to obtain y initial guess' 'equation 43 and 44 */

            if (f_y >= 0)
            {
                a1 = c1 * Tr + (v[k - 1] + v[k]);
                b1 = v[k - 1] * Tr;
            }

            else if (f_y < 0)
            {
                a1 = -c1 * Tr + (v[k - 1] + v[k]);
                b1 = -v[k] * Tr;
            }

            if (a1 <= 0)
            {
                Tao = (a1 - sqrt(pow(a1, 2) - 4 * b1 * c1)) / 2 * c1;
            }
            else if (a1 > 0)
            {
                Tao = 2 * b1 / (a1 + sqrt(pow(a1, 2) - 4 * b1 * c1));
            }

            if (f_y >= 0)
            {
                init_y = Tao + d[k - 1];
            }

            else if (f_y < 0)
            {
                init_y = Tao + d[k];
            }

            /* computing correction n to y for better approximation */

            Tr_k = d[k - 1] - init_y;
            Tr_k_1 = d[k] - init_y;

            f_ = 0;
            f__ = 0;
            for (int q = 0; q < mode; q = q + 1)
            {
                f_ = f_ + v[q] / (d[q] - init_y);
                f__ = f__ + v[q] / (pow((d[q] - init_y), 2));
            }

            f_y = f_ + p;
            f_y_ = f__;

            o_ = 0;
            for (int q = k; q < mode; q = q + 1)
            {
                o_ = o_ + v[q] / (pow((d[q] - init_y), 2));
            }
            // print(q)
            o_y_ = o_;

            u_ = 0;
            for (int q = 0; q < mode; q = q + 1)
            {
                u_ = u_ + v[q] / (pow((d[q] - init_y), 2));
            }
            // print(q)

            u_y_ = u_;
            /*'second a, b and c for final correction' */

            a2 = (Tr_k + Tr_k_1) * f_y - Tr_k * Tr_k_1 * f_y_;
            b2 = Tr_k * Tr_k_1 * f_y;

            if (f_y >= 0)
            {
                c2 = f_y - Tr_k_1 * f_y_ - u_y_ * (d[k - 1] - d[k]);
            }

            else if (f_y < 0)
            {
                c2 = f_y - Tr_k * f_y_ - o_y_ * (d[k] - d[k - 1]);
            }

            if (a2 <= 0)
            {
                n_cor = (a2 - sqrt(pow(a2, 2) - 4 * b2 * c2)) / 2 * c2;
            }

            else if (a2 > 0)
            {
                n_cor = 2 * b2 / (a2 + sqrt(pow(a2, 2) - 4 * b2 * c2));
            }

            lamda_new = init_y + n_cor;

            /*computing second correction n to y for better approximation */
            // init_y = lamda_new;

            // Tr_k = d[k - 1] - init_y;
            // Tr_k_1 = d[k] - init_y;

            // f_ = 0;
            // f__ = 0;
            // for (int q = 0; q < mode; q = q + 1)
            // {
            //     f_ = f_ + v[q] / (d[q] - init_y);
            //     f__ = f__ + v[q] / (pow((d[q] - init_y), 2));
            // }

            // f_y = f_ + p;
            // f_y_ = f__;

            // o_ = 0;
            // for (int q = k; q < mode; q = q + 1)
            // {
            //     o_ = o_ + v[q] / (pow((d[q] - init_y), 2));
            // }
            // // print(q)
            // o_y_ = o_;

            // u_ = 0;
            // for (int q = 0; q < mode; q = q + 1)
            // {
            //     u_ = u_ + v[q] / (pow((d[q] - init_y), 2));
            // }
            // // print(q)

            // u_y_ = u_;
            // /*'second a, b and c for final correction' */

            // a2 = (Tr_k + Tr_k_1) * f_y - Tr_k * Tr_k_1 * f_y_;
            // b2 = Tr_k * Tr_k_1 * f_y;

            // // c2 = f_y - Tr_k*f_y_ - o_y_*(d[k]-d[k-1]);

            // if (f_y >= 0)
            // {
            //     c2 = f_y - Tr_k_1 * f_y_ - u_y_ * (d[k - 1] - d[k]);
            // }

            // else if (f_y < 0)
            // {
            //     c2 = f_y - Tr_k * f_y_ - o_y_ * (d[k] - d[k - 1]);
            // }

            // if (a2 <= 0)
            // {
            //     n_cor = (a2 - sqrt(pow(a2, 2) - 4 * b2 * c2)) / 2 * c2;
            // }

            // else if (a2 > 0)
            // {
            //     n_cor = 2 * b2 / (a2 + sqrt(pow(a2, 2) - 4 * b2 * c2));
            // }

            // lamda_new = init_y + n_cor;

            // /*computing third correction n to y for better approximation */
            // init_y = lamda_new;

            // Tr_k = d[k - 1] - init_y;
            // Tr_k_1 = d[k] - init_y;

            // f_ = 0;
            // f__ = 0;
            // for (int q = 0; q < mode; q = q + 1)
            // {
            //     f_ = f_ + v[q] / (d[q] - init_y);
            //     f__ = f__ + v[q] / (pow((d[q] - init_y), 2));
            // }

            // f_y = f_ + p;
            // f_y_ = f__;

            // o_ = 0;
            // for (int q = k; q < mode; q = q + 1)
            // {
            //     o_ = o_ + v[q] / (pow((d[q] - init_y), 2));
            // }
            // // print(q)
            // o_y_ = o_;

            // u_ = 0;
            // for (int q = 0; q < mode; q = q + 1)
            // {
            //     u_ = u_ + v[q] / (pow((d[q] - init_y), 2));
            // }
            // // print(q)

            // u_y_ = u_;
            // /*'second a, b and c for final correction' */

            // a2 = (Tr_k + Tr_k_1) * f_y - Tr_k * Tr_k_1 * f_y_;
            // b2 = Tr_k * Tr_k_1 * f_y;

            // // c2 = f_y - Tr_k*f_y_ - o_y_*(d[k]-d[k-1]);

            // if (f_y >= 0)
            // {
            //     c2 = f_y - Tr_k_1 * f_y_ - u_y_ * (d[k - 1] - d[k]);
            // }

            // else if (f_y < 0)
            // {
            //     c2 = f_y - Tr_k * f_y_ - o_y_ * (d[k] - d[k - 1]);
            // }

            // if (a2 <= 0)
            // {
            //     n_cor = (a2 - sqrt(pow(a2, 2) - 4 * b2 * c2)) / 2 * c2;
            // }

            // else if (a2 > 0)
            // {
            //     n_cor = 2 * b2 / (a2 + sqrt(pow(a2, 2) - 4 * b2 * c2));
            // }

            // lamda_new = init_y + n_cor;

            sol[k - 1] = lamda_new;
        }

        else if (k == mode)
        {
            t = 0;
            for (int r = 0; r < mode; r = r + 1)
            {
                t = t + pow(v[r], 2);
            }

            float d6;

            d6 = d[k - 1] + t / p;

            lamda = (d6 + d[k - 1]) / 2;

            //'To check if f_y is positive or negative'

            f_ = 0;
            for (int q = 0; q < mode; q = q + 1)
            {
                f_ = f_ + v[q] / (d[q] - lamda);
            }

            f_y = f_ + p;

            //'Now, applying equation 40 and 41, choosing our initial guess'
            /*g_ = 0;
            for(int q = 0; q<n; q = q +1)
            {
                g_ = g_ + v[q]/(d[q]-lamda);
            }

            g_x = g_ + p;*/

            g_x = f_y - v[k - 1] / (d[k - 1] - lamda);
            Tr = d[k - 1] - d[k - 2];
            c1 = g_x;

            //'first a,b and c to obtain y initail guess' 'equation 43 and 44'
            if (f_y <= 0)
            {
                a1 = c1 * Tr + (v[k - 2] + v[k - 1]);
                b1 = v[k - 1] * Tr;
            }

            else if (f_y > 0)
            {
                a1 = -c1 * Tr + (v[k - 2] + v[k - 1]);
                b1 = -v[k - 1] * Tr;
            }

            if (a1 <= 0)
            {
                Tao = (a1 + sqrt(pow(a1, 2) - 4 * b1 * c1)) / 2 * c1;
            }

            else if (a1 > 0)
            {
                Tao = 2 * b1 / (a1 - sqrt(pow(a1, 2) - 4 * b1 * c1));
            }

            init_y = Tao + d[k - 1];

            //'computing correction n to y for better approximation'
            Tr_k = d[k - 1] - init_y;
            Tr_k_1 = d6 - init_y;

            f_ = 0;
            f__ = 0;
            for (int q = 0; q < mode; q = q + 1)
            {
                f_ = f_ + v[q] / (d[q] - init_y);
                f__ = f__ + v[q] / (pow((d[q] - init_y), 2));
            }

            f_y = f_ + p;
            f_y_ = f__;

            u_ = 0;
            for (int q = 0; q < mode; q = q + 1)
            {
                u_ = u_ + v[q] / (pow((d[q] - init_y), 2));
            }

            u_y_ = u_;
            //'second a, b and c for final correction'

            a2 = (Tr_k + Tr_k_1) * f_y - Tr_k * Tr_k_1 * f_y_;
            b2 = Tr_k * Tr_k_1 * f_y;

            c2 = f_y - Tr_k * u_y_ - (v[k - 1] / Tr_k_1);

            if (a2 <= 0)
            {
                n_cor = (a2 - sqrt(pow(a2, 2) - 4 * b2 * c2)) / 2 * c2;
            }

            else if (a2 > 0)
            {
                n_cor = 2 * b2 / (a2 + sqrt(pow(a2, 2) - 4 * b2 * c2));
            }

            lamda_new = init_y + n_cor;

            //'computing second correction n to y for better approximation'
            // init_y = lamda_new;

            // Tr_k = d[k - 1] - init_y;
            // Tr_k_1 = d6 - init_y;

            // f_ = 0;
            // f__ = 0;
            // for (int q = 0; q < mode; q = q + 1)
            // {
            //     f_ = f_ + v[q] / (d[q] - init_y);
            //     f__ = f__ + v[q] / (pow((d[q] - init_y), 2));
            // }

            // f_y = f_ + p;
            // f_y_ = f__;

            // u_ = 0;
            // for (int q = 0; q < mode; q = q + 1)
            // {
            //     u_ = u_ + v[q] / (pow((d[q] - init_y), 2));
            // }

            // u_y_ = u_;
            // //'second a, b and c for final correction'

            // a2 = (Tr_k + Tr_k_1) * f_y - Tr_k * Tr_k_1 * f_y_;
            // b2 = Tr_k * Tr_k_1 * f_y;

            // c2 = f_y - Tr_k * u_y_ - (v[k - 1] / Tr_k_1);

            // if (a2 <= 0)
            // {
            //     n_cor = (a2 - sqrt(pow(a2, 2) - 4 * b2 * c2)) / 2 * c2;
            // }

            // else if (a2 > 0)
            // {
            //     n_cor = 2 * b2 / (a2 + sqrt(pow(a2, 2) - 4 * b2 * c2));
            // }

            // lamda_new = init_y + n_cor;

            // //'computing third correction n to y for better approximation'
            // init_y = lamda_new;

            // Tr_k = d[k - 1] - init_y;
            // Tr_k_1 = d6 - init_y;

            // f_ = 0;
            // f__ = 0;
            // for (int q = 0; q < mode; q = q + 1)
            // {
            //     f_ = f_ + v[q] / (d[q] - init_y);
            //     f__ = f__ + v[q] / (pow((d[q] - init_y), 2));
            // }

            // f_y = f_ + p;
            // f_y_ = f__;

            // u_ = 0;
            // for (int q = 0; q < mode; q = q + 1)
            // {
            //     u_ = u_ + v[q] / (pow((d[q] - init_y), 2));
            // }

            // u_y_ = u_;
            // //'second a, b and c for final correction'

            // a2 = (Tr_k + Tr_k_1) * f_y - Tr_k * Tr_k_1 * f_y_;
            // b2 = Tr_k * Tr_k_1 * f_y;

            // c2 = f_y - Tr_k * u_y_ - (v[k - 1] / Tr_k_1);

            // if (a2 <= 0)
            // {
            //     n_cor = (a2 - sqrt(pow(a2, 2) - 4 * b2 * c2)) / 2 * c2;
            // }

            // else if (a2 > 0)
            // {
            //     n_cor = 2 * b2 / (a2 + sqrt(pow(a2, 2) - 4 * b2 * c2));
            // }

            // lamda_new = init_y + n_cor;

            sol[mode-1] = lamda_new;
        }
    }

    // gettimeofday(&tv2, NULL);

    // printf ("Total time = %f seconds\n",
    //      (double) (tv2.tv_usec - tv1.tv_usec) / 1000000 +
    //      (double) (tv2.tv_sec - tv1.tv_sec));

    // clock_t end = clock();
    // double time_spent = (double)(end-begin) / CLOCKS_PER_SEC;
    // printf("%f\n", time_spent);

    // printf("%f\n", sol[0]);
    // printf("%f\n", sol[1]);
    // printf("%f\n", sol[2]);
    // printf("%f", sol[3]);
    // float x =sol[0];
    // float s =sol[1];
    // float y =sol[2];
    // float q =sol[3];
    // printf("%f\n", y);
    // printf("%f\n", b);

    return sol;
}
