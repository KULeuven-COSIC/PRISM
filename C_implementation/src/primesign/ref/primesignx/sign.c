#include <primesign.h>
#include <curve_extras.h>
#include <tools.h>
#include <fips202.h>
#include <stdio.h>
#include <string.h>

#define RESPONSE_LENGTH TORSION_PLUS_EVEN_POWER + 16

const clock_t time_isogenies_odd = 0;
const clock_t time_sample_response = 0;
const clock_t time_change_of_basis_matrix = 0;

void
secret_sig_init(signature_t *sig)
{
    // ibz_mat_2x2_init(&(sig->mat_Baux_can_to_B_aux));
    sig->hint_aux = (int *)malloc(2 * sizeof(int));
}

void
secret_sig_finalize(signature_t *sig)
{
    // ibz_mat_2x2_finalize(&(sig->mat_Baux_can_to_B_aux));
    free(sig->hint_aux);
}

static void
ibz_vec_2_print2(char *name, const ibz_vec_2_t *vec)
{
    printf("%s", name);
    for (int i = 0; i < 2; i++) {
        ibz_printf("%Zd ", &((*vec)[i]));
    }
    ibz_printf("\n");
}

static void
ibz_vec_4_print2(char *name, const ibz_vec_4_t *vec)
{
    printf("%s", name);
    for (int i = 0; i < 4; i++) {
        ibz_printf("%Zd ", &((*vec)[i]));
    }
    ibz_printf("\n");
}



void
print_public_key(const public_key_t *pk)
{
    fp2_t j;
    ec_j_inv(&j, &pk->curve);
    fp2_print("j_EA = ", &j);
}


// compute the challenge as the hash of the message and the commitment curve and public key
void
hash_to_challenge(ibz_t *challenge_degree,
                  const unsigned char *message,
                  const public_key_t *pk,
                  size_t length)
{
    
    int found=0;
    int count=0;
    int counter_length = 4;
    char count_seq[counter_length+1];

    unsigned char *buf = malloc(FP2_ENCODED_BYTES + length + counter_length);

    {
        
        fp2_t j1;
        ec_j_inv(&j1, &pk->curve);
        fp2_encode(buf, &j1);
        memcpy(buf + FP2_ENCODED_BYTES,
               message,
               length); // TODO use defined constant
    }
    
    

    {
        ibz_t tmp;
        ibz_init(&tmp);
        ibz_pow(&tmp,&ibz_const_two,exp_size);
        digit_t digits[NWORDS_FIELD];
        while (!found) {
            
            snprintf(count_seq,counter_length+1,"%04x",count);
            
            memcpy(buf + FP2_ENCODED_BYTES + length,
               count_seq,
               counter_length);
            
            // FIXME should use SHAKE128 for smaller parameter sets?
            SHAKE256(
            (void *)digits, sizeof(digits), buf, FP2_ENCODED_BYTES + length + counter_length);

            ibz_set(challenge_degree, 1);
            ibz_copy_digit_array(challenge_degree, digits);
            ibz_mod(challenge_degree,challenge_degree,&tmp);
            found = ibz_probab_prime(challenge_degree,20);
            count++;
        }
        
    }
    free(buf);
}

int
protocols_sign(signature_t *sig,
               const public_key_t *pk,
               secret_key_t *sk,
               const unsigned char *m,
               size_t l,
               int verbose)
{
    clock_t t = tic();

    ibz_t challenge_degree;
    ec_curve_t E_aux;
    ec_basis_t Baux0; // basis of 2^TORSION_PLUS_EVEN_POWER-torsion

    quat_left_ideal_t lideal_resp;
    ibz_mat_2x2_t mat_Baux_to_Baux_can;
    // ibz_mat_2x2_t mat_sigma_phichall_BA_to_Bcomcan, mat_sigma_phichall_BA0_to_Bcom0;
    ibz_t tmp;

    ec_curve_init(&E_aux);

    ibz_init(&tmp);
    ibz_init(&challenge_degree);

    ibz_mat_2x2_init(&mat_Baux_to_Baux_can);

    quat_left_ideal_init(&lideal_resp);

    // computing the challenge
    // vec_chall is a pair of coefficients encoding the kernel of the challenge isogeny
    // as vec_chall[0]*B[0] + vec_chall[1]*B[1] where B is the canonical basis of the
    // 2^TORSION_PLUS_EVEN_TORSION torsion of EA
    hash_to_challenge(&challenge_degree, m, pk, l);

   
    // now we compute the response ideal by sampling a random ideal norm q (2^f -q)
    // computing the norm
    ibz_pow(&tmp,&ibz_const_two,exp_size);
    ibz_sub(&tmp,&tmp,&challenge_degree);
    ibz_mul(&tmp,&tmp,&challenge_degree);

    // sampling the ideal
    sampling_random_ideal_O0(&lideal_resp,&tmp,0);

    // intersection with the secret key
    quat_lideal_inter(&lideal_resp,&sk->secret_ideal,&lideal_resp,&QUATALG_PINFTY);


    // now we evaluate this isogeny on the basis of E0
    dim2id2iso_arbitrary_isogeny_evaluation(&Baux0, &E_aux, &lideal_resp);

    // notational conventions:
    // B0 = canonical basis of E0
    // Baux0 = image through isogeny of canonical basis of E0

#ifndef NDEBUG
    // testing
    assert(test_point_order_twof(&Baux0.P, &E_aux, TORSION_PLUS_EVEN_POWER));
    assert(test_point_order_twof(&Baux0.Q, &E_aux, TORSION_PLUS_EVEN_POWER));
    assert(test_point_order_twof(&Baux0.PmQ, &E_aux, TORSION_PLUS_EVEN_POWER));
#endif

    // now it only remains to format the response for the verification
    // canonical basis
    ec_basis_t B_aux_can;

    t = tic();
    
    ec_curve_to_basis_2f_to_hint(
        &B_aux_can, &E_aux, TORSION_PLUS_EVEN_POWER, sig->hint_aux);

    // compute the matrix to go from B_aux_can to Baux0
    change_of_basis_matrix_two(&mat_Baux_to_Baux_can,
                               &Baux0,                               
                               &B_aux_can,
                               &E_aux,
                               TORSION_PLUS_EVEN_POWER);

    // now we multiply by the public key matrix
    ibz_2x2_mul_mod(&mat_Baux_to_Baux_can,&mat_Baux_to_Baux_can,&sk->mat_BAcan_to_BA0_two,&TORSION_PLUS_2POWER);

    // and now we multiply by the inverse of q mod 2^f
    ibz_invmod(&tmp,&challenge_degree,&TORSION_PLUS_2POWER);
    ibz_mul(&mat_Baux_to_Baux_can[0][0],&mat_Baux_to_Baux_can[0][0],&tmp); 
    ibz_mul(&mat_Baux_to_Baux_can[0][1],&mat_Baux_to_Baux_can[0][1],&tmp);
    ibz_mul(&mat_Baux_to_Baux_can[1][0],&mat_Baux_to_Baux_can[1][0],&tmp);
    ibz_mul(&mat_Baux_to_Baux_can[1][1],&mat_Baux_to_Baux_can[1][1],&tmp);

    // computing the actual basis
    matrix_application_even_basis(&B_aux_can,&E_aux,&mat_Baux_to_Baux_can,TORSION_PLUS_EVEN_POWER);
    
    // filling the output
    copy_point(&sig->B_aux.P,&B_aux_can.P);
    copy_point(&sig->B_aux.PmQ,&B_aux_can.PmQ);
    copy_point(&sig->B_aux.Q,&B_aux_can.Q);
    
    // setting sig->E_aux
    fp2_t temp_fp2;
    fp2_copy(&temp_fp2, &E_aux.C);
    fp2_inv(&temp_fp2);
    fp2_mul(&sig->E_aux.A, &temp_fp2, &E_aux.A);
    fp2_set_one(&sig->E_aux.C);
    ec_point_init(&sig->E_aux.A24);
    sig->E_aux.is_A24_computed_and_normalized = 0;

    quat_left_ideal_finalize(&lideal_resp);

    ibz_finalize(&challenge_degree);
    ibz_finalize(&tmp);

    return 1;
}

int
protocols_verif(signature_t *sig, const public_key_t *pk, const unsigned char *m, size_t l)
{

    int verif;

    ibz_t degree_challenge;

    ibz_init(&degree_challenge);
    hash_to_challenge(&degree_challenge, m, pk, l);

    // checking that we are given A coefficients and no precomputation
    assert(fp2_is_one(&pk->curve.C) && !pk->curve.is_A24_computed_and_normalized);
    assert(fp2_is_one(&sig->E_aux.C) && !sig->E_aux.is_A24_computed_and_normalized);

    clock_t t = tic();

    // computation of the challenge
    ec_basis_t B_pk_can;
    ec_curve_t Epk;
    copy_curve(&Epk, &pk->curve);
    // ec_curve_normalize_A24(&Epk);
    ec_curve_to_basis_2f_from_hint(
        &B_pk_can, &Epk, TORSION_PLUS_EVEN_POWER, pk->hint_pk); // canonical

    ec_basis_t B_aux_can;
    ec_curve_t E_aux;
    copy_curve(&E_aux, &sig->E_aux);

    // // recovering the canonical basis
    // ec_curve_to_basis_2f_from_hint(
    //     &B_aux_can, &E_aux, TORSION_PLUS_EVEN_POWER, sig->hint_aux);

    // // TOC_clock(t,"challenge and canonical basis");

    // t = tic();

    // // applying the change matrix on the basis of E_chall
    // matrix_application_even_basis(&B_aux_can,
    //                               &E_aux,
    //                               &sig->mat_Baux_can_to_B_aux,
    //                               TORSION_PLUS_EVEN_POWER);


    // now compute the dim2 isogeny from E_pk x E_aux -> E x E'
    // of kernel B_pk_can x B_aux_can

    // first we set-up the kernel
    theta_couple_curve_t E0xEaux;
    theta_couple_point_t T1, T2, T1m2;
    theta_chain_t F;
    copy_curve(&E0xEaux.E1, &Epk);
    copy_curve(&E0xEaux.E2, &E_aux);
    copy_point(&T1.P2, &sig->B_aux.P);
    copy_point(&T2.P2, &sig->B_aux.Q);
    copy_point(&T1m2.P2, &sig->B_aux.PmQ);
    copy_point(&T1.P1, &B_pk_can.P);
    copy_point(&T2.P1, &B_pk_can.Q);
    copy_point(&T1m2.P1, &B_pk_can.PmQ);


    // multiplying by 2^(TORSION_PLUS_EVEN_POWER - exp_size + 2)
    ec_dbl_iter(&T1.P1,TORSION_PLUS_EVEN_POWER - exp_size - 2,&E0xEaux.E1,&T1.P1);
    ec_dbl_iter(&T2.P1,TORSION_PLUS_EVEN_POWER - exp_size - 2,&E0xEaux.E1,&T2.P1);
    ec_dbl_iter(&T1m2.P1,TORSION_PLUS_EVEN_POWER - exp_size - 2,&E0xEaux.E1,&T1m2.P1);
    ec_dbl_iter(&T1.P2,TORSION_PLUS_EVEN_POWER - exp_size - 2,&E0xEaux.E2,&T1.P2);
    ec_dbl_iter(&T2.P2,TORSION_PLUS_EVEN_POWER - exp_size - 2,&E0xEaux.E2,&T2.P2);
    ec_dbl_iter(&T1m2.P2,TORSION_PLUS_EVEN_POWER - exp_size - 2,&E0xEaux.E2,&T1m2.P2);


    // computing the isogeny
    int extra_info = 1;
    theta_chain_comput_strategy_faster_no_eval(
        &F,
        exp_size,
        &E0xEaux,
        &T1,
        &T2,
        &T1m2,
        strategies[TORSION_PLUS_EVEN_POWER - exp_size],
        extra_info);

    // TOC_clock(t,"response isogeny");

    // performing the degree check by evaluation and pairing
    fp2_t w0,w1,pow,invpow;
    digit_t digit_q[NWORDS_ORDER]={0};
    ec_set_zero(&T1.P2);
    ec_set_zero(&T2.P2);
    ec_set_zero(&T1m2.P2);
    copy_point(&T1.P1, &B_pk_can.P);
    copy_point(&T2.P1, &B_pk_can.Q);
    copy_point(&T1m2.P1, &B_pk_can.PmQ);

    ec_dbl_iter(&T1.P1,TORSION_PLUS_EVEN_POWER - exp_size ,&Epk,&T1.P1);
    ec_dbl_iter(&T2.P1,TORSION_PLUS_EVEN_POWER - exp_size ,&Epk,&T2.P1);
    ec_dbl_iter(&T1m2.P1,TORSION_PLUS_EVEN_POWER - exp_size ,&Epk,&T1m2.P1);

    ec_curve_normalize_A24(&Epk);
    weil(&w0,exp_size,&T1.P1,&T2.P1,&T1m2.P1,&Epk.A24);
    theta_chain_eval_special_case(&T1,&F,&T1,&E0xEaux);
    theta_chain_eval_special_case(&T2,&F,&T2,&E0xEaux);
    theta_chain_eval_special_case(&T1m2,&F,&T1m2,&E0xEaux); 
    ec_curve_normalize_A24(&F.codomain.E1);
    weil(&w1,exp_size,&T1.P1,&T2.P1,&T1m2.P1,&F.codomain.E1.A24);
    
    ibz_to_digit_array(digit_q,&degree_challenge);
    fp2_pow_vartime(&pow,&w0,digit_q,NWORDS_ORDER);
    fp2_copy(&invpow,&pow);
    fp2_inv(&invpow);
    // fp2_pow_vartime(&invpow,&w2,digit_q,&degree_challenge);

    verif = (fp2_is_equal(&pow,&w1) || fp2_is_equal(&w1,&invpow));

    ibz_finalize(&degree_challenge);

    return verif;
}
