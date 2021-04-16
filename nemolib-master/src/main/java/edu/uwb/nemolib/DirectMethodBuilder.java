package edu.uwb.nemolib;

import org.apache.commons.math3.special.Gamma;

import java.util.*;

/**
 * @Author: Zican Li
 * @Date：11/12/2020 6:22 PM
 */
public class DirectMethodBuilder {

    double EXP = 2.71828182845904523536028747135;
    long[] MSKSEG = {0x0000000000000000L,
            0xFF00000000000000L,
            0x00FF000000000000L,
            0x0000FF0000000000L,
            0x000000FF00000000L,
            0x00000000FF000000L,
            0x0000000000FF0000L,
            0x000000000000FF00L,
            0x00000000000000FFL};

    long[] DELSEG = {0x0000000000000000L,
            0x0000FFFFFFFFFFFFL,
            0xFF0000FFFFFFFFFFL,
            0xFFFF0000FFFFFFFFL,
            0xFFFFFF0000FFFFFFL,
            0xFFFFFFFF0000FFFFL,
            0xFFFFFFFFFF0000FFL,
            0xFFFFFFFFFFFF0000L};

    long MSKBIT[] = {0x0000000000000000L,
            0x8080808080808080L,
            0x4040404040404040L,
            0x2020202020202020L,
            0x1010101010101010L,
            0x0808080808080808L,
            0x0404040404040404L,
            0x0202020202020202L,
            0x0101010101010101L};

    long DELBIT[] = {0x0000000000000000L,
            0x3F3F3F3F3F3F3F3FL,
            0x9F9F9F9F9F9F9F9FL,
            0xCFCFCFCFCFCFCFCFL,
            0xE7E7E7E7E7E7E7E7L,
            0xF3F3F3F3F3F3F3F3L,
            0xF9F9F9F9F9F9F9F9L,
            0xFCFCFCFCFCFCFCFCL};

    boolean directed;
    long n;
    long[] r;
    long[] c;
    double b;
    long f;
    double[] logfac_r;
    double[] logfac_c;
    double[] logfac_f;
    double[] logfac_cache;
    long sqsum_r;
    long sqsum_c;
    long maxdeg;
    Random rand;

    List<Short> odg;
    List<Short> idg;

    public DirectMethodBuilder(long n, long[] r, long[] c, Random rand) {
        this.directed = true;
        this.n = n;
        this.r = r;
        this.c = c;
        long tmp = 0L;
        this.rand = rand;
        f = 0;
        sqsum_r = 0;
        sqsum_c = 0;
        logfac_r = new double[(int) n];
        logfac_c = new double[(int) n];
        maxdeg = 0;
        for (int i = 0; i != n; i++) {
            tmp += (c[i]) * (r[i]);
            logfac_r[i] = Gamma.logGamma((double) (r[i] + 1));
            logfac_c[i] = Gamma.logGamma((double) (c[i] + 1));
            if (r[i] > maxdeg) {
                maxdeg = r[i];
            }
            if (c[i] > maxdeg) {
                maxdeg = c[i];
            }
            f += r[i];
            sqsum_r += r[i] * r[i] - r[i];
            sqsum_c += c[i] * c[i] - c[i];
        }
        logfac_f = new double[64];
        for (int i = 0; i != 64; ++i) {
            logfac_f[i] = Gamma.logGamma((double) (f - i + 1));
        }
        logfac_cache = new double[(int) (maxdeg + 1)];
        for (int i = 0; i <= maxdeg; ++i) {
            logfac_cache[i] = Gamma.logGamma(i + 1);
        }
        this.b = (double) (tmp) / (double) (f);

    }

    public DirectMethodBuilder(long n, long[] r, Random rand) {
        this.directed = false;
        this.n = n;
        this.r = r;
        this.b = 0;
        this.rand = rand;
        f = 0;
        sqsum_r = 0;
        logfac_r = new double[(int) n];
        for (int i = 0; i != n; ++i) {
            logfac_r[i] = Gamma.logGamma((double) (r[i] + 1));
            f += r[i];
            sqsum_r += r[i] * r[i] - r[i];
        }
    }

    long SET(long g, long row, long col) {
        g |= (1L << (63 - (row * 8 + col)));
        return g;
    }

    void noautoms(long g, short k, boolean directed) {

        long gr = 0L;
        Set<Long> already = new HashSet<>();
        short[] odeg = new short[k];
        short[] ideg = new short[k];
        for (int i = 0; i != k; ++i) {
            odeg[i] = 0;
            ideg[i] = 0;
        }
        //Convert graph format
        for (int i = 0; i != k; ++i) {
            for (int j = 0; j != k; ++j) {
                if ((g & (1L << i + j * k)) > 0) {
                    gr = SET(gr, j, i);
                    ++ideg[i];
                    ++odeg[j];
                }
            }
        }

        long tmp1;
        long tmp2;
        short[] c = new short[k + 1];
        short[] o = new short[k + 1];
        short j = k;
        short s = 0;
        short q;
        short t1;
        short t2;
        for (int i = 0; i != k + 1; ++i) {
            c[i] = 0;
            o[i] = 1;
        }
        already.add(gr);
        for (int i = 0; i != k; ++i) {
            idg.add(ideg[i]);
            odg.add(odeg[i]);
        }
        while (j != 0) {
            q = (short) (c[j] + o[j]);
            if (q < 0) {
                o[j] = (short) -o[j];
                --j;
            } else if (q == j) {
                ++s;
                o[j] = (short) -o[j];
                --j;
            } else if (q > -1) {
                t1 = (short) (j - q + s);
                t2 = (short) (j - c[j] + s);
                if (t1 > t2) {
                    t1 ^= t2;
                    t2 ^= t1;
                    t1 ^= t2;
                }
                tmp1 = ideg[t1 - 1];
                ideg[t1 - 1] = ideg[t2 - 1];
                ideg[t2 - 1] = (short) tmp1;
                tmp1 = odeg[t1 - 1];
                odeg[t1 - 1] = odeg[t2 - 1];
                odeg[t2 - 1] = (short) tmp1;
                tmp1 = gr & MSKSEG[t1];
                tmp2 = gr & MSKSEG[t2];
                gr &= DELSEG[t1];
                gr |= (tmp1 >>> 8);
                gr |= (tmp2 << 8);
                tmp1 = gr & MSKBIT[t1];
                tmp2 = gr & MSKBIT[t2];
                gr &= DELBIT[t1];
                gr |= (tmp1 >>> 1);
                gr |= (tmp2 << 1);
                if (!already.contains(gr)) {
                    already.add(gr);
                    for (int i = 0; i != k; ++i) {
                        idg.add(ideg[i]);
                        odg.add(odeg[i]);
                    }
                }
                c[j] = q;
                j = k;
                s = 0;
            }
        }
    }

    long getcanonical(long g, short k) {
        long gr = g;
        long tmp1;
        long tmp2;
        short[] c = new short[k + 1];
        short[] o = new short[k + 1];
        short j = k;
        short s = 0;
        short q;
        short t1;
        short t2;
        long canon = g;

        for (int i = 0; i != k + 1; ++i) {
            c[i] = 0;
            o[i] = 1;
        }

        while (j != 0) {
            q = (short) (c[j] + o[j]);
            if (q < 0) {
                o[j] = (short) -o[j];
                --j;
            } else if (q == j) {
                ++s;
                o[j] = (short) -o[j];
                --j;
            } else if (q > -1) {
                t1 = (short) (j - q + s);
                t2 = (short) (j - c[j] + s);
                if (t1 > t2) {
                    t1 ^= t2;
                    t2 ^= t1;
                    t1 ^= t2;
                }
                tmp1 = gr & MSKSEG[t1];
                tmp2 = gr & MSKSEG[t2];
                gr &= DELSEG[t1];
                gr |= (tmp1 >> 8);
                gr |= (tmp2 << 8);
                tmp1 = gr & MSKBIT[t1];
                tmp2 = gr & MSKBIT[t2];
                gr &= DELBIT[t1];
                gr |= (tmp1 >> 1);
                gr |= (tmp2 << 1);
                // compare to max
                if (gr > canon) {
                    canon = gr;
                }
                //continue algorithm
                c[j] = q;
                j = k;
                s = 0;
            }
        }
        return canon;
    }

    List<Double> calc_motif_undir(long[] positions, short[] degrees, short k, List<Double> result) {
        long fn = this.f;
        long sqn = this.sqsum_r;
        for (int i = 0; i != k; ++i) {
            if (r[(int) positions[i]] < degrees[i]) {
                //cout << "GA";
                return result;
            }
            fn -= r[(int) positions[i]];
            sqn = sqn - (r[(int) positions[i]] * r[(int) positions[i]] - r[(int) positions[i]]);
        }
        double olda = (double) (this.sqsum_r) / (double) (2 * this.f);
        double newa = (double) (sqn) / (double) (2 * fn);
        double oldb = 0;
        long be = 0;
        for (int i = 0; i != k; ++i)
            for (int j = i + 1; j < k; ++j)
                be += r[(int) positions[i]] * r[(int) positions[i]];
        double newb = (double) (be) / (double) (fn);
        double changelog = 0;
        for (int i = 0; i != k; ++i)
            changelog += logfac_r[(int) positions[i]]
                    - Gamma.logGamma(r[(int) positions[i]] - degrees[i] + 1);
        double x = (((double) (fn) * Math.log((double) (fn) / EXP) - (double) (f) * Math.log((double) (f) / EXP)) / 2
                + olda * olda + olda - newa * newa - newa - newb
                + changelog);
        result.add(x);
        return result;
    }

    List<Double> calc_motif_dir(long[] positions, short[] degrees_r, short[] degrees_c, short k, List<Double> result) {
        long delta_f = 0;
        long sqn_r = this.sqsum_r;
        long sqn_c = this.sqsum_c;
        for (int i = 0; i != k; ++i) {
            if (r[(int) positions[i]] < degrees_r[i] || c[(int) positions[i]] < degrees_c[i]) {
                return result;
            }
            delta_f += degrees_r[i];
            long rold = r[(int) positions[i]];
            long rnew = r[(int) positions[i]] - degrees_r[i];
            sqn_r = sqn_r - (rold * rold - rold) + (rnew * rnew - rnew);
            long cold = c[(int) positions[i]];
            long cnew = c[(int) positions[i]] - degrees_c[i];
            sqn_c = sqn_c - (cold * cold - cold) + (cnew * cnew - cnew);
        }

        double olda = 0.5 * ((double) (this.sqsum_r) / (double) (this.f))
                * ((double) (this.sqsum_c) / (double) (this.f));
        double newa = 0.5 * ((double) (sqn_r) / (double) (this.f - delta_f))
                * ((double) (sqn_c) / (double) (this.f - delta_f));

        double oldb = this.b;
        long be = 0L;
        for (int i = 0; i != k; ++i)
            for (int j = 0; j != k; ++j) {
                be += (r[(int) positions[i]] - degrees_r[i]) * (c[(int) positions[j]] - degrees_c[i]);
                if (i == j)
                    be -= (r[(int) positions[i]]) * (c[(int) positions[j]]);
            }
        double newb = (oldb * (double) (this.f) + (double) (be)) / (double) (this.f - delta_f);

        double changelog = 0;
        for (int i = 0; i != k; ++i)
            changelog += logfac_r[(int) positions[i]]
                    - logfac_cache[(int) (r[(int) positions[i]] - degrees_r[i])]
                    + logfac_c[(int) positions[i]]
                    - logfac_cache[(int) (c[(int) positions[i]] - degrees_c[i])];
        double x = (logfac_f[(int) delta_f] - logfac_f[0] + olda + oldb - newa - newb + changelog);
        result.add(x);
        return result;
    }

    public double motif_sample_log(long motif, short k, long num_samples) {
        List<Double> results = new ArrayList<>();
        List<Double> result_buffer = new ArrayList<>();
        long BUF_SIZE = 1000000; //maximum size of results vector, should be at least 1e6\
        odg = new ArrayList<>();
        idg = new ArrayList<>();
        noautoms(motif, k, this.directed);

        int num_graphs = odg.size() / k;
        short[] odeg = new short[k];
        short[] ideg = new short[k];
        long[] pos = new long[k];

        // Figure out the suited vertices
        List<Long> candidate_indices = new ArrayList<>();
        for (int i = 0; i != n; ++i) {
            if (directed) {
                for (int j = 0; j != k; ++j) {
                    if ((r[i] >= odg.get(j)) && (c[i] >= idg.get(j))) {
                        candidate_indices.add((long) i);
                        break;
                    }
                }
            } else {
                for (int j = 0; j != k; ++j) {
                    if ((r[i] >= odg.get(j))) {
                        candidate_indices.add((long) i);
                        break;
                    }
                }
            }
        }

        long new_n = candidate_indices.size();

        //cout << "For motif ID " << motif << " we have " << new_n
        //	 << " candidates out of " << n << " vertices." << endl;
        boolean[] taken = new boolean[(int) new_n];
        for (int i = 0; i != new_n; ++i) {
            taken[i] = false;
        }
        double maximum = 0.0;
        long cand;
        if (new_n >= k) {
            for (int i = 0; i != num_samples; ++i) {
                for (int j = 0; j != k; ++j) {
                    do {
                        cand = Math.floorMod((rand.nextInt() * 10) ^ 8,new_n);
                    } while (taken[(int) cand]);
                    taken[(int) cand] = true;
                    pos[j] = cand;
                }

                for (int j = 0; j != k; ++j) {
                    taken[(int) pos[j]] = false;
                    pos[j] = candidate_indices.get((int) pos[j]);
                }

                for (int j = 0; j != num_graphs; ++j) {
                    if (directed) {
                        for (int l = 0; l != k; ++l) {
                            odeg[l] = odg.get(j * k + l);
                            ideg[l] = idg.get(j * k + l);
                        }
                        results = calc_motif_dir(pos, odeg, ideg, k, results);
                    } else {
                        for (int l = 0; l != k; ++l) {
                            odeg[l] = odg.get(j * k + l);
                        }
                        results = calc_motif_undir(pos, odeg, k, results);
                    }

                    if (results.size() > BUF_SIZE) {

                        double sum = 0.0;
                        long l = 0;
                        Collections.sort(results);
                        if (result_buffer.size() == 0)
                            maximum = results.get(results.size() - 1);
                        while (l != BUF_SIZE) {
                            sum += Math.exp(results.get((int) l) - maximum);
                            ++l;
                        }
                        results = results.subList((int) BUF_SIZE, results.size());
                        result_buffer.add(sum);
                    }
                }
            }
        }
        //apply correctional factor due to vertex selection
        double p = 1.0;

        for (int i = 0; i != k; ++i) {
            p *= (double) (new_n - i) / (double) (n - i);
        }
        //cout << "correctional factor " << p << endl;

        //Calculate results
        Collections.sort(results);
        double ret = 0.0;
        if (results.size() > 0) {
            //cout << "Number of actual samples: "
            //	 << (result_buffer.size() * BUF_SIZE + results.size()) << endl;
            if (result_buffer.size() == 0)
                maximum = results.get(results.size() - 1);
            //double maximum = results[results.size()-1];
            double sum = 0;
            for (long i = 0L; i != results.size(); ++i) {
                sum += Math.exp(results.get((int) i) - maximum);
            }
            for (long i = 0L; i != result_buffer.size(); ++i) {
                sum += result_buffer.get((int) i);
            }

            ret = maximum + Math.log(sum) - Math.log((double) (num_samples) / (p));
        } else {
            ret = Math.sqrt(-1.0);
        }
        return ret;
    }

    public static long[] getDegrees(Graph graph) {
        int size = graph.getSize();
        long[] degrees = new long[size];
        for (int i = 1; i <= size; i++) {
            degrees[i-1] = graph.getAdjacencyList(i).size();
        }
        return degrees;
    }

}
