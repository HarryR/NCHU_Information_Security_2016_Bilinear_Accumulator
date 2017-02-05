#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include <pbc.h>
#include <pbc_test.h>
#include <pbc_poly.h>
#include "./pbc/misc/darray.c"
#include "./pbc/arith/fp.c"
#include "./pbc/arith/poly.c"

#define SIZEOFELEMENT 1000

// Query by these values
#define MEMBERSHIP 94
#define NONMEMBERSHIP 87

/* Function declaration */
void poly_extgcd(element_ptr, element_ptr, element_ptr, element_ptr, element_ptr );
void element_generate(pairing_t, field_t, element_t *);
void element_set_x(element_t, int);
void accumulate(pairing_t, field_t, element_t, element_t, element_t *, element_t );
int membership_verify(pairing_t, field_t, element_t, element_t, element_t, int );
int nonmembership_verify(pairing_t, field_t, element_t, element_t, element_t, int);
void ask_server_for_Wx_and_Xj(pairing_t, field_t, element_t, element_t, element_t, element_t, int);
void ask_server_for_Ax_and_Bx(pairing_t, field_t, element_t, element_t, element_t, element_t, element_t, element_t);

void server_emulator(pairing_t, field_t, element_t, element_t, element_t *);
void respond_Wx_and_Xj(pairing_t, field_t, element_t, element_t, element_t, element_t, int );
void respond_Ax_and_Bx(pairing_t, field_t, element_t, element_t, element_t, element_t, element_t,  element_t);

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

        // Clear initilized variables
        element_clear(temp_u);
        element_clear(temp_v);
}

void element_generate(pairing_t pairing, field_t fx, element_t *poly_x) {
        for (int i = 0; i < SIZEOFELEMENT; i++) {
                element_init(poly_x[i], fx);
                poly_alloc(poly_x[i], 1);
                element_ptr constcoeff;
                constcoeff = (element_ptr)poly_coeff(poly_x[i], 0);
                // element_random(constcoeff); // Randomly generate
                element_set_si(constcoeff, i); // Set ith element as i
        }
}

void element_set_x(element_t X, int x) {
        poly_alloc(X, 1);
        element_ptr constcoeff;
        constcoeff = (element_ptr)poly_coeff(X, 0);
        element_set_si(constcoeff, x);
}

void accumulate(pairing_t pairing, field_t fx, element_t g, element_t secret_key, element_t *poly_x, element_t acc) {
        element_t xi_plus_s;
        element_t x0_mul_to_xn;

        element_init(xi_plus_s, fx);
        element_init(x0_mul_to_xn, fx);

        poly_set0(xi_plus_s);
        poly_set1(x0_mul_to_xn);

        // Calculate x0_mul_to_xn = (x0+s)(x1+s)...(xn+s)
        for (int i = 0; i < SIZEOFELEMENT; i++) {
                poly_add(xi_plus_s, poly_x[i], secret_key);
                // element_printf("add = %B\n", xi_plus_s);
                poly_mul(x0_mul_to_xn, xi_plus_s, x0_mul_to_xn);
                // element_printf("mul = %B\n", x0_mul_to_xn);
        }

        // Calculate g^(x0+s)(x1+s)...(xn+s)
        element_pow_zn(acc, g, x0_mul_to_xn);

        // Clear initilized variables
        element_clear(xi_plus_s);
        element_clear(x0_mul_to_xn);
}

int membership_verify(pairing_t pairing, field_t fx, element_t g, element_t secret_key, element_t acc, int j) {
        element_t temp5, temp6;
        element_t temp_add, temp_mul;
        element_t h, public_key;

        element_init_G1(public_key, pairing);
        element_init_G1(h, pairing);

        element_init_GT(temp5, pairing);
        element_init_GT(temp6, pairing);

        // Calculate e(acc(X), g)
        element_pairing(temp5, acc, g);
        // element_printf("\ntemp5 = %B\n", temp5);

        // Ask for Wx and (Xj+s) and assign then to public_key and temp_add
        ask_server_for_Wx_and_Xj(pairing, fx, g, secret_key, public_key, temp_add, j);

        // Calculate g^(xj+s) and assign it to h
        element_pow_zn(h, g, temp_add);

        // Calculate e(Wx, g^(xj+s))
        element_pairing(temp6, public_key, h);
        // element_printf("\ntemp6 = %B\n", temp6);

        // Check if e(acc(X), g) equals e(Wx, g^(xj+s))
        int verified;
        if (!element_cmp(temp5, temp6)) {
                verified = 1;
        } else {
                verified = 0;

        }

        // Clear initilized variables
        element_clear(h);
        element_clear(public_key);

        element_clear(temp_add);

        element_clear(temp5);
        element_clear(temp6);

        return verified;
}

int nonmembership_verify(pairing_t pairing, field_t fx, element_t g, element_t secret_key, element_t acc, int x) {
        element_t Ax;
        element_t Bx;
        element_t sig;
        element_t temp1, temp2, temp3, temp4;

        element_t X;
        element_t h;

        element_t x_plus_s;
        element_t eea_gcd;

        element_init(x_plus_s, fx);
        element_init(eea_gcd, fx);

        element_init(X, fx);
        element_init_G1(h, pairing);

        element_init_G1(Ax, pairing);
        element_init_G1(Bx, pairing);
        element_init_G1(sig, pairing);

        element_init_GT(temp1, pairing);
        element_init_GT(temp2, pairing);
        element_init_GT(temp3, pairing);
        element_init_GT(temp4, pairing);

        // Set X
        element_set_x(X, x);
        // element_printf("X = %B\n", X);

        // Ask for GCD and Ax and Bx and assign then to eea_gcd and Ax and Bx
        ask_server_for_Ax_and_Bx(pairing, fx, g, secret_key, eea_gcd, X, Ax, Bx);

        // Calculate e(acc(X), Ax)
        element_pairing(temp3, acc, Ax);
        // element_printf("temp3 = %B\n", temp3);

        // Calculate (x+s)
        poly_add(x_plus_s, X, secret_key);
        // element_printf("X+S = %B\n", x_plus_s);

        // Calculate g^(x+s)
        element_pow_zn(h, g, x_plus_s);
        // element_printf("h = %B\n", h);

        // Calculate e(g^(x+s), Bx),  where h = g^(x+s)
        element_pairing(temp4, h, Bx);
        // element_printf("temp4 = %B\n", temp4);

        // Calculate e(acc(X), Ax) * e(g^(x+s), Bx)
        element_mul(temp2, temp3, temp4);

        // Calculate e(g, g)
        int eea_gcd_degree = poly_coeff_count(eea_gcd) - 1;
        element_pow_zn(sig, g, element_item(eea_gcd, eea_gcd_degree));
        element_pairing(temp1, sig, g);
        // element_printf("\ntemp1 = %B\n", temp1);
        // element_printf("\ntemp2 = %B\n", temp2);

        int verified;
        if (!element_cmp(temp1, temp2)) {
                verified = 1;
        } else {
                verified = 0;

        }

        // Clear initilized variables
        element_clear(x_plus_s);

        element_clear(X);
        element_clear(eea_gcd);

        element_clear(h);

        element_clear(Ax);
        element_clear(Bx);

        element_clear(sig);
        element_clear(temp1);
        element_clear(temp2);
        element_clear(temp3);
        element_clear(temp4);

        return verified;
}

void ask_server_for_Wx_and_Xj(pairing_t pairing, field_t fx, element_t g, element_t secret_key, element_t public_key, element_t temp_add, int j) {
        respond_Wx_and_Xj(pairing, fx, g, secret_key, public_key, temp_add, j);
}

void ask_server_for_Ax_and_Bx(pairing_t pairing, field_t fx, element_t g, element_t secret_key, element_t eea_gcd, element_t X, element_t Ax, element_t Bx) {
        respond_Ax_and_Bx( pairing,  fx,  g,  secret_key,  eea_gcd,  X, Ax,Bx);
}

/* Start of Server logic */
element_t poly_x_server[SIZEOFELEMENT]; // Array of element_t saved on server
void server_emulator(pairing_t pairing, field_t fx, element_t g, element_t secret_key, element_t *poly_x) {
        for (int i = 0; i < SIZEOFELEMENT; i++) {
                element_init(poly_x_server[i], fx);
                poly_alloc(poly_x_server[i], 1);
                element_set(poly_x_server[i], poly_x[i]);
        }
}

void respond_Wx_and_Xj(pairing_t pairing, field_t fx, element_t g, element_t secret_key, element_t public_key, element_t temp_add, int j) {
        element_t temp_mul;

        element_init(temp_add, fx);
        element_init(temp_mul, fx);

        // Calculate Wx = g^(x1+s)(x2+s)...(xn+s) => where (xj+s) is exlcluded
        poly_set0(temp_add);
        poly_set1(temp_mul);
        for (int i = 0; i < SIZEOFELEMENT; i++) {
                // (x1+s)(x2+s)...(xn+s) => where (xj+s) is exlcluded
                if (i != j) {
                        poly_add(temp_add, poly_x_server[i], secret_key);
                        poly_mul(temp_mul, temp_add, temp_mul);
                }
        }
        element_pow_zn(public_key, g, temp_mul);

        // Calculate (xj+s)
        poly_add(temp_add, poly_x_server[j], secret_key);// (xj+s)

        // Clear initilized variables
        element_clear(temp_mul);
}

void respond_Ax_and_Bx(pairing_t pairing, field_t fx, element_t g, element_t secret_key, element_t eea_gcd, element_t X, element_t Ax, element_t Bx) {
        element_t eea_x, eea_y;

        element_t x_plus_s;
        element_t x0_mul_to_xn;

        element_init(eea_x, fx);
        element_init(eea_y, fx);

        element_init(x_plus_s, fx);
        element_init(x0_mul_to_xn, fx);

        // Calculate (x0+s)(x1+s)...(xn+s)
        poly_set0(x_plus_s);
        poly_set1(x0_mul_to_xn);
        for (int i = 0; i < SIZEOFELEMENT; i++) {
                poly_add(x_plus_s, poly_x_server[i], secret_key);
                // element_printf("add = %B\n", xi_plus_s);
                poly_mul(x0_mul_to_xn, x_plus_s, x0_mul_to_xn);
                // element_printf("mul = %B\n", x0_mul_to_xn);
        }
        // element_printf("(x0+s)(x1+s)...(xn+s) = %B\n", x0_mul_to_xn);

        // Calculate (x+s)
        poly_add(x_plus_s, X, secret_key);
        // element_printf("X+S = %B\n", x_plus_s);

        // Calculate Alpha(s) and Beta(s) and assign them to eea_x and eea_y
        poly_extgcd(eea_gcd, eea_x, eea_y, x0_mul_to_xn, x_plus_s);
        // element_printf("gcd = %B\n", eea_gcd);
        // element_printf("x = %B\n", eea_x);
        // element_printf("y = %B\n", eea_y);

        // Calculate Ax = g^Alpha(s)
        element_pow_zn(Ax, g, eea_x);
        // element_printf("k = %B\n", Ax);

        // Calculate Bx = g^Beta(s)
        element_pow_zn(Bx, g, eea_y);
        // element_printf("Bx = %B\n", Bx);

        // Clear initilized variables
        element_clear(eea_x);
        element_clear(eea_y);
}
/* End of Server logic */

int main(int argc, char **argv) {

        pairing_t pairing;
        field_t fx;

        element_t g;
        element_t secret_key;

        element_t poly_x[SIZEOFELEMENT];
        element_t acc;

        pbc_demo_pairing_init(pairing, argc, argv);

        field_init_fp(pairing->Zr, pairing->r);
        field_init_poly(fx, pairing->Zr);

        element_init(secret_key, fx);
        element_init_G1(g, pairing);

        element_init_G1(acc, pairing);


        printf("\nGenerating elements...\n");
        element_generate(pairing, fx, poly_x);
        for (int i = 0; i < SIZEOFELEMENT; i++) {
                element_printf("%B ", poly_x[i]);
        }
        printf("\n");

        printf("\nGenerating g...\n");
        element_random(g);
        element_printf("%B\n", g);

        printf("\nGenerating s...\n");
        int polyX_degree = 2;
        do {
                poly_random_monic(secret_key, polyX_degree);
        } while (!poly_is_irred(secret_key));
        element_printf("%B\n", secret_key);


        // Pass g, secret_key and elements (x0...xn) to server
        printf("\nPassing elements to server...\n");
        server_emulator(pairing, fx, g, secret_key, poly_x);
        for (int i = 0; i < SIZEOFELEMENT; i++) {
                element_printf("%B ", poly_x_server[i]);
        }
        printf("\n");


        printf("Accumulating g^(x0+s)(x1+s)...(xn+s)...\n");
        // Accumulate g^(x0+s)(x1+s)...(xn+s)
        accumulate(pairing, fx, g, secret_key, poly_x, acc);
        // element_printf("acc = %B\n", acc);

        printf("\nVerifying by membership witness... ");
        int j = MEMBERSHIP;
        if (membership_verify(pairing, fx, g, secret_key, acc, j)) {
                printf("signature verifies\n");
                printf("===>\t%d is a member\n", j);
        } else {
                printf("signature does not verify\n");
                printf("===>\t%d is not a member\n", j);
        }

        printf("\nVerifying by non-membership witness... ");
        int x = NONMEMBERSHIP;
        if (nonmembership_verify(pairing, fx, g, secret_key, acc, x)) {
                printf("signature verifies\n");
                printf("===>\t%d is not a member\n", x);
        } else {
                printf("signature does not verify\n");
                printf("===>\t%d is a member\n", x);
        }


        // Clear initilized variables
        for (int i = 0; i < SIZEOFELEMENT; i++) {
                element_clear(poly_x[i]);
        }

        for (int i = 0; i < SIZEOFELEMENT; i++) {
                element_clear(poly_x_server[i]);
        }

        element_clear(acc);

        element_clear(secret_key);
        element_clear(g);

        field_clear(fx);
        pairing_clear(pairing);

        return 0;
}
