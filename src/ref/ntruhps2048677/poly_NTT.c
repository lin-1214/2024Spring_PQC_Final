#include "poly_NTT.h"
#include "params.h"

#define Q_mask 0x03FF
#define sz 1536
#define newQ 1419778561
#define newZ 4360860
#define mod_newQ 422821617
// #define DEBUG
#define uint16_t uint64_t

typedef struct {
    uint64_t low;
    uint64_t high;
} uint128_t;

uint128_t multiply_uint64(uint64_t a, uint64_t b) {
    uint128_t result;
    __uint128_t product = (__uint128_t)a * (__uint128_t)b;
    result.low = (uint64_t)product;
    result.high = (uint64_t)(product >> 64);
    return result;
}

uint16_t mod_inverse(int a, int m){
    for (int x = 1; x < m; x++)
        if (((a % m) * (x % m)) % m == 1)
            return x;
}


void good_thomas_permutation(uint16_t **r1, uint16_t **r2, const poly *a, const poly *b) {
    // printf("[+] Good Thomas permutation start.\n");
    uint16_t pad_a[sz]={0};
    uint16_t pad_b[sz]={0};
    for (size_t i = 0; i < NTRU_N; i++) {
        pad_a[i] = a->coeffs[i];
        pad_b[i] = b->coeffs[i];
    }
    for (size_t i = NTRU_N; i < sz; i++) {
        pad_a[i] = 0;
        pad_b[i] = 0;
    }
    uint16_t idx = 0;
    for (size_t i = 0; i < 3; i++) {
        for (size_t j = 0; j < 512; j++) {
            idx = (i * 512 + j * 3) % sz;     
            r1[i][j] = pad_a[idx];
            r2[i][j] = pad_b[idx];
        }
    }
}

void inv_good_thomas_permutation(uint16_t *r, uint16_t **c) {
    
    uint16_t idx = 0;
    for (size_t i = 0; i < 512; i++) {
        r[idx] = c[0][i];
        idx = (idx + 3) % sz;
    }

    idx = 512;
    for (size_t i = 0; i < 512; i++) {
        r[idx] = c[1][i];
        idx = (idx + 3) % sz;
    }

    idx = 1024;
    for (size_t i = 0; i < 512; i++) {
        r[idx] = c[2][i];
        idx = (idx + 3) % sz;
    }
    // for (size_t i = 0; i < sz; i++) {
    //     printf("%lld ", r[i]);
    // }
}

void ntt_512(uint16_t *p, uint16_t zetas[], uint16_t brv[]) {
    //  for (size_t i=0; i<512; i++) {
    //          printf("%lld ", p[i]);
    //  }
    //  printf("\n");
     uint16_t window_sz = 256;
     while (window_sz >= 1) {
        uint16_t start = 0;
        uint16_t brv_idx = 0;
        while (start < 512) {
            uint16_t zeta = zetas[brv[brv_idx]];
            // printf("zeta: %lld\n", zeta);
            uint16_t i = start;
            while (i < start + window_sz) {
                uint128_t res = multiply_uint64(p[i + window_sz], zeta);
                uint16_t b = res.low % newQ;
                bool overflow = (res.high != 0);
                while (overflow) {
                    uint128_t tmp = multiply_uint64((res.high%newQ), mod_newQ);
                    if (tmp.high == 0) {
                        overflow = false;
                    }
                    else {
                        res.high = tmp.high;
                        b = (b + tmp.low%newQ) % newQ;
                    }
                }
                // printf("b: %lld i: %lld\n", b, i);
                // if (i == 384) {
                //     printf("%lld, %lld\n", res.high, res.low);
                //     printf("%lld, %lld\n", p[i+window_sz], zeta);
                //     printf("%lld\n", b);
                //     printf("%lld\n", i+window_sz);
                // }

                // if (i+window_sz == 448) {
                //     printf("%lld, %lld\n", res.high, res.low);
                //     printf("%lld, %lld\n", p[i+window_sz], zeta);
                //     printf("%lld\n", b);
                //     printf("%lld, %lld\n", i, window_sz);
                // }
                uint16_t tmp = p[i]%newQ;
                while (tmp < b%newQ) {
                    tmp += newQ;
                }
                p[i+window_sz] = (tmp - b%newQ) % newQ;

                // if (i>400) exit(0);
                p[i] = (p[i]%newQ + b%newQ) % newQ;
                // printf("p[i]: %lld, p[i+window_sz]: %lld\n", p[i], p[i+window_sz]);
                i += 1;
            }
            start = window_sz + i;
            brv_idx += 1;
        }
        window_sz /= 2;     // window_sz = window_sz / 2
     }
}

void intt_512(uint16_t *p, uint16_t *invzetas, uint16_t *brv) {
    uint16_t window_sz = 1;
    uint16_t inv_2 = mod_inverse(2, newQ); 
    while (window_sz <= 256) {
        uint16_t start = 0;
        uint16_t brv_idx = 0;
        while (start < 512) {
            uint16_t zeta = invzetas[brv[brv_idx]];
            uint16_t i = start;
            while (i < start + window_sz) {
                // deal with a
                 // uint16_t a = (p[i] + p[i+window_sz]) * inv_2;
                uint16_t tmp = (p[i]%newQ + p[i+window_sz]%newQ)%newQ;
                uint128_t res = multiply_uint64(tmp, inv_2);
                uint16_t a = res.low % newQ;
                bool overflow = (res.high != 0);
                while (overflow) {
                    uint128_t tmp = multiply_uint64((res.high%newQ), mod_newQ);
                    if (tmp.high == 0) {
                        overflow = false;
                    }
                    else {
                        res.high = tmp.high;
                        a = (a + tmp.low%newQ) % newQ;
                    }
                }

                // uint16_t b = (p[i] - p[i+window_sz]) * inv_2;
                tmp = p[i]%newQ;
                while (tmp < p[i+window_sz]%newQ) {
                    tmp += newQ;
                }

                tmp = (tmp - p[i+window_sz]%newQ) % newQ;
                res = multiply_uint64(tmp, inv_2);
                uint16_t b = res.low % newQ;
                overflow = (res.high != 0);
                while (overflow) {
                    uint128_t tmp = multiply_uint64((res.high%newQ), mod_newQ);
                    if (tmp.high == 0) {
                        overflow = false;
                    }
                    else {
                        res.high = tmp.high;
                        b = (b + tmp.low%newQ) % newQ;
                    }
                }

                p[i] = a % newQ;
                res = multiply_uint64(b, zeta);
                p[i+window_sz] = res.low % newQ;
                overflow = (res.high != 0);
                while (overflow) {
                    uint128_t tmp = multiply_uint64((res.high%newQ), mod_newQ);
                    if (tmp.high == 0) {
                        overflow = false;
                    }
                    else {
                        res.high = tmp.high;
                        p[i+window_sz] = (p[i+window_sz] + tmp.low%newQ) % newQ;
                    }
                }
                i += 1;
            }
            start = window_sz + i;
            brv_idx += 1;
        }
        window_sz = window_sz << 1;      // window_sz = window_sz * 2
    }
}

void conv(uint16_t *r, uint16_t *a, uint16_t *b) {
    for (size_t i = 0; i < 3; i++) {
        for (size_t j = 0; j < 3; j++) {
            uint16_t idx = 0;
            idx = (i + j) % 3;

            uint128_t res = multiply_uint64(a[i], b[j]);
            r[idx] += res.low % newQ;
            bool overflow = (res.high != 0);
            while (overflow) {
                uint128_t tmp = multiply_uint64((res.high%newQ), mod_newQ);
                if (tmp.high == 0) {
                    overflow = false;
                }
                else {
                    res.high = tmp.high;
                    r[idx] = (r[idx] + tmp.low%newQ) % newQ;
                }
            }

            // r[idx] += (a[i] * b[j]) % newQ;
            if (r[idx] >= newQ) {
                r[idx] %= newQ;
            }
        }
    }
}

void final_stage(poly *r, uint16_t *p) {
    for (size_t i = 0; i < NTRU_N; i++) {
        r->coeffs[i] = 0;
    }

    for (size_t i = 0; i < sz; i++) {
        uint16_t idx = i % NTRU_N;
        // printf("%lld\n ", idx);
        uint16_t tmp = (r->coeffs[idx])%NTRU_Q;
        r->coeffs[idx] = (tmp + p[i]%NTRU_Q) % NTRU_Q;
    }

    for (size_t i = 0; i < NTRU_N; i++) {
        r->coeffs[i] = r->coeffs[i] % NTRU_Q;
    }
}

void poly_NTT(poly *r, const poly *a, const poly *b) {
    uint16_t zetas[] = {1, 4360860, 585893566, 392536624, 66488721, 826132640, 439732608, 1205515279, 1132394312, 319036243, 143596738, 578317142, 
         335065015, 836145506, 1276692374, 1070111631, 745047005, 286996363, 540466509, 1059363934, 1272963317, 1318240305, 739430276, 
         537903795, 499443525, 507600816, 759785221, 1047722775, 1410553693, 1067659455, 920937097, 1245350916, 342752416, 704288034, 
         1038330015, 90886382, 244276882, 268877342, 962803782, 993026782, 1220989152, 1290757884, 96879348, 866015315, 328110608, 
         110898446, 524882935, 975672559, 340714355, 819387434, 512922246, 490435476, 100264424, 1270841958, 225956480, 1339796214, 
         922006767, 534115670, 140213260, 43251974, 1061065912, 758427489, 645152625, 1277565510, 847867989, 837040388, 760546705, 
         910784836, 191535558, 1186466458, 294895240, 1189552308, 593525082, 755930544, 99364234, 1355999723, 1197062098, 1124056456, 
         625267293, 1363608748, 25295906, 789563704, 797016290, 330160399, 36658650, 433806083, 688849662, 938994154, 364611193, 
         1257749275, 1263199715, 1356088414, 366549205, 714985962, 700892196, 554620565, 43072619, 1217007723, 1308311413, 1218237973, 
         913144394, 758403871, 1294196098, 805179228, 1071885687, 781516276, 1176619008, 986222807, 570938430, 666359199, 61731293, 
         152997892, 169022146, 1036105288, 1134587036, 272810548, 778725501, 1339377400, 368016173, 434862576, 725483197, 1172750534, 
         921556969, 1161980692, 3392168, 96917421, 782946458, 912877055, 568650912, 484783305, 549434885, 251507427, 1221063354, 
         1165779696, 951725860, 1111410448, 193686775, 445901990, 212972849, 891947673, 928656277, 558964406, 239975017, 1391778496, 
         852125783, 616562470, 48787864, 387680068, 142484998, 979606557, 770880072, 925852877, 1037229177, 1181579004, 508002532, 
         580276707, 76532256, 527344451, 356636598, 1221470709, 21146624, 9642568, 307447343, 1010800094, 1014933921, 1127356124, 
         559746038, 1359152137, 717371975, 212469807, 754088298, 849349373, 952818078, 764125773, 1099542243, 1220893169, 137229438, 
         283756619, 686907180, 445834560, 54345176, 847030679, 1115812680, 885203014, 547679213, 877488980, 958594307, 378893451, 
         1339457207, 549422748, 1275334925, 481886300, 305846363, 412341731, 815756550, 106634278, 1344810433, 1219457146, 485202229, 
         163341518, 909017536, 1221465227, 251147641, 1099997299, 805973051, 1402721749, 1249411031, 1297195646, 994188015, 845402445, 
         927875001, 992253446, 785042567, 1216395077, 248349094, 25058674, 1272595153, 156353195, 1357591621, 524211088, 148135482, 
         272558081, 133868656, 758049302, 1207279687, 268864572, 645985461, 1365205432, 346621002, 1257370192, 757553339, 704819910, 
         556934701, 1299956552, 1196529895, 179450550, 918885337, 1251728738, 1089367907, 446035459, 143386179, 956725369, 1150645472, 
         1126061305, 494698795, 1257081030, 525916187, 465836665, 918281880, 509909251, 1291641831, 167200653, 1411542, 796984185, 
         882827638, 569876909};

    uint16_t invzetas[] = {1, 849901652, 536950923, 622794376, 1418367019, 1252577908, 128136730, 909869310, 501496681, 953941896, 893862374, 162697531, 
                925079766, 293717256, 269133089, 463053192, 1276392382, 973743102, 330410654, 168049823, 500893224, 1240328011, 223248666, 119822009, 
                862843860, 714958651, 662225222, 162408369, 1073157559, 54573129, 773793100, 1150913989, 212498874, 661729259, 1285909905, 1147220480, 
                1271643079, 895567473, 62186940, 1263425366, 147183408, 1394719887, 1171429467, 203383484, 634735994, 427525115, 491903560, 574376116, 
                425590546, 122582915, 170367530, 17056812, 613805510, 319781262, 1168630920, 198313334, 510761025, 1256437043, 934576332, 200321415, 
                74968128, 1313144283, 604022011, 1007436830, 1113932198, 937892261, 144443636, 870355813, 80321354, 1040885110, 461184254, 542289581, 
                872099348, 534575547, 303965881, 572747882, 1365433385, 973944001, 732871381, 1136021942, 1282549123, 198885392, 320236318, 655652788, 
                466960483, 570429188, 665690263, 1207308754, 702406586, 60626424, 860032523, 292422437, 404844640, 408978467, 1112331218, 1410135993, 
                1398631937, 198307852, 1063141963, 892434110, 1343246305, 839501854, 911776029, 238199557, 382549384, 493925684, 648898489, 440172004, 
                1277293563, 1032098493, 1370990697, 803216091, 567652778, 28000065, 1179803544, 860814155, 491122284, 527830888, 1206805712, 
                973876571, 1226091786, 308368113, 468052701, 253998865, 198715207, 1168271134, 870343676, 934995256, 851127649, 506901506, 636832103, 
                1322861140, 1416386393, 257797869, 498221592, 247028027, 694295364, 984915985, 1051762388, 80401161, 641053060, 1146968013, 285191525, 
                383673273, 1250756415, 1266780669, 1358047268, 753419362, 848840131, 433555754, 243159553, 638262285, 347892874, 614599333, 125582463, 
                661374690, 506634167, 201540588, 111467148, 202770838, 1376705942, 865157996, 718886365, 704792599, 1053229356, 63690147, 156578846, 
                162029286, 1055167368, 480784407, 730928899, 985972478, 1383119911, 1089618162, 622762271, 630214857, 1394482655, 56169813, 794511268, 
                295722105, 222716463, 63778838, 1320414327, 663848017, 826253479, 230226253, 1124883321, 233312103, 1228243003, 508993725, 659231856, 
                582738173, 571910572, 142213051, 774625936, 661351072, 358712649, 1376526587, 1279565301, 885662891, 497771794, 79982347, 1193822081, 
                148936603, 1319514137, 929343085, 906856315, 600391127, 1079064206, 444106002, 894895626, 1308880115, 1091667953, 553763246, 
                1322899213, 129020677, 198789409, 426751779, 456974779, 1150901219, 1175501679, 1328892179, 381448546, 715490527, 1077026145, 
                174427645, 498841464, 352119106, 9224868, 372055786, 659993340, 912177745, 920335036, 881874766, 680348285, 101538256, 146815244, 
                360414627, 879312052, 1132782198, 674731556, 349666930, 143086187, 583633055, 1084713546, 841461419, 1276181823, 1100742318, 
                287384249, 214263282, 980045953, 593645921, 1353289840, 1027241937, 833884995, 1415417701};

    uint16_t brv[] = {0, 128, 64, 192, 32, 160, 96, 224, 16, 144, 80, 208, 48, 176, 112, 240, 8, 136, 72, 200, 40, 168, 104, 232, 24, 152, 88, 216, 56, 184, 120, 
        248, 4, 132, 68, 196, 36, 164, 100, 228, 20, 148, 84, 212, 52, 180, 116, 244, 12, 140, 76, 204, 44, 172, 108, 236, 28, 156, 92, 220, 60, 
        188, 124, 252, 2, 130, 66, 194, 34, 162, 98, 226, 18, 146, 82, 210, 50, 178, 114, 242, 10, 138, 74, 202, 42, 170, 106, 234, 26, 154, 90, 
        218, 58, 186, 122, 250, 6, 134, 70, 198, 38, 166, 102, 230, 22, 150, 86, 214, 54, 182, 118, 246, 14, 142, 78, 206, 46, 174, 110, 238, 30, 
        158, 94, 222, 62, 190, 126, 254, 1, 129, 65, 193, 33, 161, 97, 225, 17, 145, 81, 209, 49, 177, 113, 241, 9, 137, 73, 201, 41, 169, 105, 
        233, 25, 153, 89, 217, 57, 185, 121, 249, 5, 133, 69, 197, 37, 165, 101, 229, 21, 149, 85, 213, 53, 181, 117, 245, 13, 141, 77, 205, 45, 
        173, 109, 237, 29, 157, 93, 221, 61, 189, 125, 253, 3, 131, 67, 195, 35, 163, 99, 227, 19, 147, 83, 211, 51, 179, 115, 243, 11, 139, 75, 
        203, 43, 171, 107, 235, 27, 155, 91, 219, 59, 187, 123, 251, 7, 135, 71, 199, 39, 167, 103, 231, 23, 151, 87, 215, 55, 183, 119, 247, 15, 
        143, 79, 207, 47, 175, 111, 239, 31, 159, 95, 223, 63, 191, 127, 255};

    uint16_t **reorder_a = (uint16_t **)malloc(3 * sizeof(uint16_t *));
    uint16_t **reorder_b = (uint16_t **)malloc(3 * sizeof(uint16_t *));
    uint16_t **ntt_conv_result = (uint16_t **)malloc(3 * sizeof(uint16_t *));
    for (size_t i = 0; i < 3; i++) {
        reorder_a[i] = (uint16_t *)malloc(512 * sizeof(uint16_t));
        reorder_b[i] = (uint16_t *)malloc(512 * sizeof(uint16_t));
        ntt_conv_result[i] = (uint16_t *)malloc(512 * sizeof(uint16_t));
    }

    uint16_t *tem1 = (uint16_t *)malloc(3 * sizeof(uint16_t));
    uint16_t *tem2 = (uint16_t *)malloc(3 * sizeof(uint16_t));
    uint16_t *conv_result = (uint16_t *)malloc(3 * sizeof(uint16_t));
    uint16_t *ntt_order = (uint16_t *)malloc(sz * sizeof(uint16_t));

    // printf("[+] Initialization done.\n");

    good_thomas_permutation(reorder_a, reorder_b, a, b);
    // printf("[+] Good Thomas permutation success.\n");

    for (size_t i = 0; i < 3; i++) {
        ntt_512(reorder_a[i], zetas, brv);
        ntt_512(reorder_b[i], zetas, brv);
    }

    // printf("[+] NTT512 success.\n");

    for (size_t i = 0; i < 512; i++) {
        for (size_t j = 0; j < 3; j++) {
            tem1[j] = reorder_a[j][i];
            tem2[j] = reorder_b[j][i];
            conv_result[j] = 0;
        }
        conv(conv_result, tem1, tem2);
        for (size_t j = 0; j < 3; j++) {
            ntt_conv_result[j][i] = conv_result[j];
        }
    }

    // free memory
    free(tem1);
    free(tem2);
    free(*reorder_a);
    free(*reorder_b);
    free(reorder_a);
    free(reorder_b);
    free(conv_result);

    // printf("[+] Convolution success.\n");

    for (size_t i = 0; i < 3; i++) {
        intt_512(ntt_conv_result[i], invzetas, brv);
    }

    // printf("[+] INTT512 success.\n");

    inv_good_thomas_permutation(ntt_order, ntt_conv_result);
    # ifdef DEBUG
    for (size_t i = 0; i<sz; i++) {
        printf("%lld ", ntt_order[i]);
    }
    printf("\n");
    # endif

    free(*ntt_conv_result);
    free(ntt_conv_result);
    // printf("[+] Inverse Good Thomas permutation success.\n");

    final_stage(r, ntt_order);
    // printf("[+] Final stage success.\n");  
    free(ntt_order);
}
