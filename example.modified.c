#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include <pbc.h>
#include <pbc_test.h>
#include <pbc_poly.h>
#include "./pbc/misc/darray.c"
#include "./pbc/arith/fp.c"
#include "./pbc/arith/poly.c"

#define SIZEOFELEMENT 3

void poly_extgcd(element_ptr d, element_ptr s, element_ptr t, element_ptr f, element_ptr g) {
        element_t a, b, q, r, u[8], v[8], temp_u, temp_v;
        element_init(a, d->field);
        element_init(b, d->field);
        element_init(q, d->field);
        element_init(r, d->field);

        element_init(temp_u, s->field);
        element_init(temp_v, t->field);

        int i;
        for (i = 0; i < 8; i++) {
                element_init(u[i], s->field);
                element_init(v[i], t->field);
        }

        element_set(a, f);
        element_set(b, g);
        element_set0(u[1]);
        element_set1(u[0]);
        element_set0(v[0]);
        element_set1(v[1]);

        for (i=1;; i++) {
                poly_div(q, r, a, b);

                if (element_is0(r)) {
                        break;
                }

                poly_mul(temp_u, q, u[i]);
                poly_sub(u[i+1], u[i-1], temp_u);

                poly_mul(temp_v, q, v[i]);
                poly_sub(v[i+1], v[i-1], temp_v);

                element_set(a, b);
                element_set(b, r);
        }
        element_set(d, b);
        element_set(s, u[i]);
        element_set(t, v[i]);

        element_clear(a);
        element_clear(b);
        element_clear(q);
        element_clear(r);

        for (i = 0; i < 8; i++) {
                element_clear(u[i]);
                element_clear(v[i]);
        }

        element_clear(temp_u);
        element_clear(temp_v);
}

int main(int argc, char **argv) {
        int polyX_degree = 2;

        pairing_t pairing;
        field_t fx;
        element_t g;
        element_t acc, k;
        element_t h, public_key;
        element_t sig;
        element_t temp1, temp2, temp3, temp4;
        element_t temp5, temp6;
        element_t secret_key, poly_x[SIZEOFELEMENT];
        element_t temp_add, temp_mul;
        element_t randomX;
        element_t eea_gcd, eea_x, eea_y;

        pbc_demo_pairing_init(pairing, argc, argv);

        field_init_fp(pairing->Zr, pairing->r);
        field_init_poly(fx, pairing->Zr);

        element_init(secret_key, fx);
        element_init(temp_add, fx);
        element_init(temp_mul, fx);
        element_init(randomX, fx);
        element_init(eea_gcd, fx);
        element_init(eea_x, fx);
        element_init(eea_y, fx);

        element_init_G1(g, pairing);
        element_init_G1(public_key, pairing);
        element_init_G1(k, pairing);
        element_init_G1(h, pairing);
        element_init_G1(sig, pairing);
        element_init_G1(acc, pairing);

        element_init_GT(temp1, pairing);
        element_init_GT(temp2, pairing);
        element_init_GT(temp3, pairing);
        element_init_GT(temp4, pairing);

        element_init_GT(temp5, pairing);
        element_init_GT(temp6, pairing);

        do {
                poly_random_monic(secret_key, polyX_degree);
        } while (!poly_is_irred(secret_key));

        element_printf("secret key = %B\n", secret_key);

        // Generating random element set
        for (int i = 0; i < SIZEOFELEMENT; i++) {
                element_init(poly_x[i], fx);
                poly_alloc(poly_x[i], 1);
                element_ptr constcoeff;
                constcoeff = (element_ptr)poly_coeff(poly_x[i], 0);
                element_random(constcoeff);
        }

        poly_set0(temp_add);
        poly_set1(temp_mul);

        // Accumulate
        printf("\nAccumulating...\n");
        for (int i = 0; i < SIZEOFELEMENT; i++) {
                poly_add(temp_add, poly_x[i], secret_key);
                element_printf("add = %B\n", temp_add);
                poly_mul(temp_mul, temp_add, temp_mul);
                element_printf("mul = %B\n", temp_mul);
        }

        element_init(randomX, fx);
        poly_alloc(randomX, 1);
        element_ptr constcoeff;
        constcoeff = (element_ptr)poly_coeff(randomX, 0);
        element_random(constcoeff);
        element_printf("ranX = %B\n", randomX);

        poly_add(temp_add, randomX, secret_key);
        element_printf("ranX+S = %B\n", temp_add);

        poly_extgcd(eea_gcd, eea_x, eea_y, temp_mul, temp_add);

        element_printf("gcd = %B\n", eea_gcd);
        element_printf("x = %B\n", eea_x);
        element_printf("y = %B\n", eea_y);

        // Generate non-membership witness
        printf("\nGenerating non-membership witness\n");
        element_random(g);

        element_pow_zn(acc, g, temp_mul); // g^(x1+s)(x2+s)...(xn+s) => (x1+s)(x2+s)...(xn+s) = temp_mul
        element_printf("\nacc = %B\n", acc);

        element_pow_zn(k, g, eea_x); // A = g^a(s) => a(s) = eea_x
        element_printf("\nk = %B\n", k);

        element_pairing(temp4, acc, k); // e(acc(X), Ax)
        element_printf("\ntemp4 = %B\n", temp4);

        element_pow_zn(h, g, temp_add); // g^(x+s) => (x+s) = temp_add
        element_printf("\nh = %B\n", h);

        element_pow_zn(public_key, g, eea_y); //public = Bx = g^b(s) => b(s) = eea_y
        element_printf("\npublic_key = %B\n", public_key);

        element_pairing(temp3, h, public_key); // e(g^x, g^s, Bx)
        element_printf("\ntemp3 = %B\n", temp3);
        element_mul(temp2, temp3, temp4); // temp4 * temp3

        int eea_gcd_degree = poly_coeff_count(eea_gcd) - 1;
        element_pow_zn(sig, g, element_item(eea_gcd, eea_gcd_degree));

        element_pairing(temp1, sig, g); // e(g, g)

        element_printf("\ntemp1 = %B\n", temp1);
        element_printf("\ntemp2 = %B\n", temp2);

        if (!element_cmp(temp1, temp2)) {
                printf("\nsignature verifies\n");
        } else {
                printf("\nsignature does not verify\n");
        }


        // Generate membership witness
        printf("\nGenerating membership witness\n");
        element_pairing(temp5, acc, g); // e(acc(X), g)
        element_printf("\ntemp5 = %B\n", temp5);

        poly_set0(temp_add);
        poly_set1(temp_mul);

        int n = 1;

        for (int i = 0; i < 3; i++) {
                if (i != n) { // W0
                        poly_add(temp_add, poly_x[i], secret_key);
                        poly_mul(temp_mul, temp_add, temp_mul); // Wx = (x1+s)(x2+s)...(xn+s) = temp_mul
                }
        }
        element_pow_zn(public_key, g, temp_mul); // g^Wx = g^temp_mull

        poly_add(temp_add, poly_x[n], secret_key);
        element_pow_zn(h, g, temp_add); // g^(x+s) => (x+s) = temp_add

        element_pairing(temp6, public_key, h); // e(Wx, g^s * g^x) // h = g^(xn+s)
        element_printf("\ntemp6 = %B\n", temp6);

        if (!element_cmp(temp5, temp6)) {
                printf("\nsignature verifies\n");
        } else {
                printf("\nsignature does not verify\n");
        }

        field_clear(fx);
        element_clear(secret_key);
        element_clear(temp_add);
        element_clear(temp_mul);
        element_clear(randomX);
        element_clear(eea_gcd);
        element_clear(eea_x);
        element_clear(eea_y);

        element_clear(sig);
        element_clear(public_key);
        element_clear(g);
        element_clear(h);
        element_clear(temp1);
        element_clear(temp2);

        element_clear(k);
        element_clear(acc);
        element_clear(temp3);
        element_clear(temp4);

        element_clear(temp5);
        element_clear(temp6);

        for (int i = 0; i < 3; i++) {
                element_clear(poly_x[i]);
        }

        pairing_clear(pairing);

        return 0;
}
