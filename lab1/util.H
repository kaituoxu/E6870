
//  $Id: util.H,v 1.26 2009/09/18 06:16:48 stanchen Exp stanchen $

/** * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * **
*   @file util.H
*   @brief I/O routines, GmmSet, and Graph classes.
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef _UTIL_H
#define _UTIL_H


#include <cassert>
#include <cfloat>
#include <cmath>
#include <algorithm>
#include <fstream>
#include <iostream>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>
#include <boost/format.hpp>
#include <boost/numeric/ublas/matrix.hpp>


using namespace std;

#ifndef SWIG
using boost::format;
using boost::str;
using boost::numeric::ublas::matrix;
#endif


/** * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * **
*   @name Math stuff.
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

//@{

/** This value can be used to represent the logprob of a "zero" prob.
*   Theoretically, log 0 is negative infinity which we can't store,
*   so we can use this very large negative value instead.
**/
const double g_zeroLogProb = -FLT_MAX / 2.0;

/** Add the log probs held in @p logProbList, returning answer as log prob.
*   That is, let's say we have a list of probability values, the logs of
*   which are stored in @p logProbList.  Then, this routine returns the
*   log of the sum of those probability values.
*   Logarithms are base <i>e</i>.
**/
double add_log_probs(const vector<double>& logProbList);

/** Does in-place real FFT.
*   For inputs <tt>vals[i]</tt>, i = 0, ..., N-1 with sample period T,
*   on return the real and imaginary parts of the FFT value for frequency
*   i/NT are held in the outputs <tt>vals[2*i]</tt> and <tt>vals[2*i+1]</tt>.
**/
void real_fft(vector<double>& vals);

/** Sets @p vec to be equal to the @p rowIdx-th row of @p mat.
*   Rows are numbered starting from 0.
**/
void copy_matrix_row_to_vector(const matrix<double>& mat, int rowIdx,
    vector<double>& vec);

/** Sets the @p rowIdx-th row of @p mat to @p vec; sizes must match.
*   Rows are numbered starting from 0.
**/
void copy_vector_to_matrix_row(const vector<double>& vec,
    matrix<double>& mat, int rowIdx);

//@}


/** * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * **
*   @name Command-line parsing and parameter lookup routines.
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

//@{

/** Type of object used for holding program parameters.
*   Declaration needed for hack to get default arguments to work.
**/
typedef map<string, string> ParamsType;

#ifndef SWIG
/** Given cmd line arguments @p argv, parses flags of the
*   form <tt>--&lt;flag&gt; &lt;val&gt;</tt> and places in flag-to-value
*   map @p params.
*   Expects same @p argv value as passed to <tt>main()</tt>.
*   Existing values in @p params are not erased (unless overriden
*   in @p argv).
**/
void process_cmd_line(const char** argv, map<string, string>& params);
#endif

/** Like process_cmd_line(), but expects arguments as string vector. **/
void process_cmd_line(const vector<string>& argList,
    map<string, string>& params);

/** Like process_cmd_line(), but expects space-separated arguments
*   in single string.
**/
void process_cmd_line(const string& argStr, map<string, string>& params);

/** Returns value of boolean parameter @p name from parameter map @p params.
*   If not present, returns @p defaultVal.
**/
bool get_bool_param(const map<string, string>& params, const string& name,
    bool defaultVal = false);

/** Like get_bool_param(), but for integer parameters. **/
int get_int_param(const map<string, string>& params, const string& name,
    int defaultVal = 0);

/** Like get_bool_param(), but for floating-point parameters. **/
double get_float_param(const map<string, string>& params, const string& name,
    double defaultVal = 0.0);

/** Like get_bool_param(), but for string parameters. **/
string get_string_param(const map<string, string>& params, const string& name,
    const string& defaultVal = string());

/** Like get_string_param(), but throws exception if parameter absent. **/
string get_required_string_param(const map<string, string>& params,
    const string& name);

//@}


/** * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * **
*   @name Vector/matrix I/O routines.
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

//@{

/** Splits @p inStr into space-separated tokens; places in @p outList. **/
void split_string(const string& inStr, vector<string>& outList);

/** Reads a list of strings, one to a line, from file @p fileName and
*   places in @p strList.
**/
void read_string_list(const string& fileName, vector<string>& strList);

/** Reads matrix of floating-point numbers from stream @p inStrm in Matlab
*   text format and places in @p mat.  Expects optional matrix header, and
*   then one row per line.  If argument @p name is provided, checks
*   name associated with matrix matches and throws exception if doesn't.
*   Returns name given in matrix header, or empty string if none provided.
**/
string read_float_matrix(istream& inStrm, matrix<double>& mat,
    const string& name = string());

/** Like read_float_matrix(), but for float vectors. **/
string read_float_vector(istream& inStrm, vector<double>& vec,
    const string& name = string());

/** Like read_float_matrix(), but for integer matrices. **/
string read_int_matrix(istream& inStrm, matrix<int>& mat,
    const string& name = string());

/** Like read_float_matrix(), but for integer vectors. **/
string read_int_vector(istream& inStrm, vector<int>& vec,
    const string& name = string());

/** Like read_float_matrix(), but reads from file @p fileName instead
*   of stream.
**/
void read_float_matrix(const string& fileName, matrix<double>& mat);

/** Like read_float_vector(), but reads from file @p fileName instead
*   of stream.
**/
void read_float_vector(const string& fileName, vector<double>& vec);

/** Like read_int_matrix(), but reads from file @p fileName instead
*   of stream.
**/
void read_int_matrix(const string& fileName, matrix<int>& mat);

/** Like read_int_vector(), but reads from file @p fileName instead
*   of stream.
**/
void read_int_vector(const string& fileName, vector<int>& vec);

/** Writes floating-point matrix @p mat to stream @p outStrm in Matlab text
*   format.  If the argument @p name is provided, this name will
*   be written in the matrix header (and is the name
*   that will be assigned to the matrix if loaded in octave).
**/
void write_float_matrix(ostream& outStrm, const matrix<double>& mat,
    const string& name = string());

/** Like write_float_matrix(), but for float vectors. **/
void write_float_vector(ostream& outStrm, const vector<double>& vec,
    const string& name = string());

/** Like write_float_matrix(), but for integer matrices. **/
void write_int_matrix(ostream& outStrm, const matrix<int>& mat,
    const string& name = string());

/** Like write_float_matrix(), but for integer vectors. **/
void write_int_vector(ostream& outStrm, const vector<int>& vec,
    const string& name = string());

/** Like write_float_matrix(), but writes to file @p fileName instead
*   of stream.
**/
void write_float_matrix(const string& fileName, const matrix<double>& mat);

/** Like write_float_vector(), but writes to file @p fileName instead
*   of stream.
**/
void write_float_vector(const string& fileName, const vector<double>& vec);

/** Like write_int_matrix(), but writes to file @p fileName instead
*   of stream.
**/
void write_int_matrix(const string& fileName, const matrix<int>& mat);

/** Like write_int_vector(), but writes to file @p fileName instead
*   of stream.
**/
void write_int_vector(const string& fileName, const vector<int>& vec);

//@}


/** * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * **
*   Class holding set of diagonal covariance GMM's.
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
class GmmSet
    {
    public:
        /** Ctor; load from file @p fileName if present. **/
        GmmSet(const string& fileName = string())
            { if (!fileName.empty()) read(fileName); }

        /** Read GMM parameters from file @p fileName. **/
        void read(const string& fileName);

        /** Write GMM parameters to file @p fileName. **/
        void write(const string& fileName) const;

        /** Make empty object, where @p gmmGaussCounts holds the number
        *   of Gaussians in each GMM, and @p dimCnt is the dimension of
        *   the Gaussians.  Component weights are initialized to be
        *   uniform; and means are variances for each dimension are set
        *   to 0 and 1, respectively.
        **/
        void init(const vector<int>& gmmGaussCounts, int dimCnt);

        /** Clear out object. **/
        void clear();

        /** Returns whether object is empty. **/
        bool empty() const { return m_gmmMap.empty(); }

        /** Recompute normalization constants. **/
        void compute_norms();

        /** Given input feature vectors, compute log prob of each GMM
        *   for each frame.
        **/
        void calc_gmm_probs(const matrix<double>& feats,
            matrix<double>& logProbs) const;

        /** Compute log prob of each Gaussian for GMM @p gmmIdx
        *   for frame @p frmIdx in features @p feats.
        *   Places result in @p logProbs.  Returns total log prob of GMM.
        **/
        double calc_gaussian_probs(const matrix<double>& feats, int frmIdx,
            int gmmIdx, vector<double>& logProbs) const;

        /** Returns number of leaves/GMM's in model. **/
        int get_gmm_count() const { return m_gmmMap.size(); }

        /** Returns total number of individual Gaussian in model. **/
        int get_gaussian_count() const { return m_gaussParams.size1(); }

        /** Return number of dimensions in expected feature vector. **/
        int get_dim_count() const
            { assert(!(m_gaussParams.size2() & 1));
                return m_gaussParams.size2() / 2; }

        /** Returns index of first Gaussian for GMM @p gmmIdx. **/
        int get_min_gaussian_index(int gmmIdx) const
            {
            assert((gmmIdx >= 0) && (gmmIdx < (int) m_gmmMap.size()));
            return m_gmmMap[gmmIdx];
            }

        /** Returns one past index of last Gaussian for GMM @p gmmIdx. **/
        int get_max_gaussian_index(int gmmIdx) const
            {
            assert((gmmIdx >= 0) && (gmmIdx < (int) m_gmmMap.size()));
            return (gmmIdx + 1 < (int) m_gmmMap.size()) ?
                m_gmmMap[gmmIdx + 1] : m_gaussParams.size1();
            }

        /** Returns mixture weight of Gaussian with index @p gaussIdx. **/
        double get_gaussian_weight(int gaussIdx) const
            {
            assert((gaussIdx >= 0) &&
                (gaussIdx < (int) m_gaussWeights.size()));
            return m_gaussWeights[gaussIdx];
            }

        /** Set mixture weight of Gaussian with index @p gaussIdx. **/
        void set_gaussian_weight(int gaussIdx, double wgt)
            {
            assert((gaussIdx >= 0) &&
                (gaussIdx < (int) m_gaussWeights.size()));
            m_gaussWeights[gaussIdx] = wgt;
            }

        /** Returns mean for dimension @p dimIdx for Gaussian @p gaussIdx. **/
        double get_gaussian_mean(int gaussIdx, int dimIdx) const
            {
            assert((gaussIdx >= 0) &&
                (gaussIdx < (int) m_gaussParams.size1()) && (dimIdx >= 0) &&
                (2 * dimIdx < (int) m_gaussParams.size2()));
            return m_gaussParams(gaussIdx, 2 * dimIdx);
            }

        /** Returns var for dimension @p dimIdx for Gaussian @p gaussIdx. **/
        double get_gaussian_var(int gaussIdx, int dimIdx) const
            {
            assert((gaussIdx >= 0) &&
                (gaussIdx < (int) m_gaussParams.size1()) && (dimIdx >= 0) &&
                (2 * dimIdx + 1 < (int) m_gaussParams.size2()));
            return m_gaussParams(gaussIdx, 2 * dimIdx + 1);
            }

        /** Sets mean for dimension @p dimIdx for Gaussian @p gaussIdx. **/
        void set_gaussian_mean(int gaussIdx, int dimIdx, double val)
            {
            assert((gaussIdx >= 0) &&
                (gaussIdx < (int) m_gaussParams.size1()) && (dimIdx >= 0) &&
                (2 * dimIdx < (int) m_gaussParams.size2()));
            m_gaussParams(gaussIdx, 2 * dimIdx) = val;
            }

        /** Sets var for dimension @p dimIdx for Gaussian @p gaussIdx. **/
        void set_gaussian_var(int gaussIdx, int dimIdx, double val)
            {
            assert((gaussIdx >= 0) &&
                (gaussIdx < (int) m_gaussParams.size1()) && (dimIdx >= 0) &&
                (2 * dimIdx + 1 < (int) m_gaussParams.size2()));
            m_gaussParams(gaussIdx, 2 * dimIdx + 1) = val;
            }

        /** Returns log norm constant + log weight for Gaussian @p gaussIdx. **/
        double get_gaussian_norm(int gaussIdx) const
            {
            assert((gaussIdx >= 0) && (gaussIdx < (int) m_logNorms.size()));
            return m_logNorms[gaussIdx];
            }

    private:
        /** For each GMM, index of first Gaussian in m_gaussParams.
        *   Assumes that Gaussians for each GMM are contiguous.
        **/
        vector<int> m_gmmMap;

        /** For each Gaussian, its mixture weight.
        *   We assume each Gaussian belongs to single mixture.
        **/
        vector<double> m_gaussWeights;

        /** For each Gaussian, alternating mean + var for each dim. **/
        matrix<double> m_gaussParams;

        /** Log norm constant + log weight for each Gaussian. **/
        vector<double> m_logNorms;
    };


/** * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * **
*   Graph symbol table and I/O routines.
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
class SymbolTable
    {
    public:
        /** Ctor; load from file @p fileName if present. **/
        SymbolTable(const string& fileName = string())
            { if (!fileName.empty()) read(fileName); }

        /** Read symbols from file @p fileName. **/
        void read(const string& fileName);

        /** Clear out object. **/
        void clear();

        /** Returns number of symbols in table. **/
        int size() const { return m_strToIdxMap.size(); }

        /** Returns whether symbol table is empty. **/
        bool empty() const { return m_strToIdxMap.empty(); }

        /** Map from a string to its index.
        *   If not in table, returns -1.
        **/
        int get_index(const string& theStr) const
            {
            map<string, int>::const_iterator lookup =
                m_strToIdxMap.find(theStr);
            return (lookup != m_strToIdxMap.end()) ? lookup->second : -1;
            }

        /** Map from an index to its string.
        *   If not in table, returns empty string.
        **/
        string get_str(int theIdx) const
            {
            map<int, string>::const_iterator lookup =
                m_idxToStrMap.find(theIdx);
            return (lookup != m_idxToStrMap.end()) ? lookup->second :
                string();
            }

    private:
        /** Map from strings to integer indices. **/
        map<string, int> m_strToIdxMap;

        /** Map from integer indices to strings. **/
        map<int, string> m_idxToStrMap;
    };


/** * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * **
*   Arc class.
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
class Arc
    {
    public:
        /** Ctor, no args. **/
        Arc() : m_dst(-1), m_gmmIdx(-1), m_wordIdx(-1), m_cost(0.0) { }

        /** Ctor, fill in all vals. **/
        Arc(int dst, int gmmIdx, int wordIdx, double cost) :
            m_dst(dst), m_gmmIdx(gmmIdx), m_wordIdx(wordIdx), m_cost(cost)
            { }

        /** Sets all values. **/
        void assign(int dst, int gmmIdx, int wordIdx, double cost)
            { m_dst = dst; m_gmmIdx = gmmIdx; m_wordIdx = wordIdx;
                m_cost = cost; }

        /** Returns dest state index.
        *   To find src state, see get_src_state() method for Graph class.
        **/
        int get_dst_state() const { return m_dst; }

        /** Returns assoc GMM index, or -1 if not present. **/
        int get_gmm() const { return  m_gmmIdx; }

        /** Returns assoc word index, or -1 if not present. **/
        int get_word() const { return  m_wordIdx; }

        /** Returns assoc cost, or neg log prob. **/
        double get_cost() const { return m_cost; }

    private:
        /** Destination state. **/
        int m_dst;

        /** GMM index, or -1 if not present. **/
        int m_gmmIdx;

        /** Word index, or -1 if not present. **/
        int m_wordIdx;

        /** Cost. **/
        float m_cost;
    };


/** * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * **
*   Graph class and I/O routines.
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
class Graph
    {
    public:
        /** Ctor; load from file @p fileName if present,
        *   and load sym table from file @p symFile if present.
        **/
        Graph(const string& fileName = string(),
            const string& symFile = string()) :
            m_symTable(new SymbolTable), m_start(-1)
            {
            if (!symFile.empty())
                read_word_sym_table(symFile);
            if (!fileName.empty())
                read(fileName);
            }

        /** Read graph from file @p fileName.
        *   Optionally, read symbol table from file @p symFile.
        **/
        void read(const string& fileName, const string& symFile = string());

        /** Read graph from stream @p inStrm.
        *   If argument @p name is provided, checks name associated with
        *   graph matches.  Returns name given in graph header,
        *   or empty string if none provided.
        **/
        string read(istream& inStrm, const string& name = string());

        /** Read word symbol table from file @p symFile.
        *   Pass in the empty string to load an empty symbol table.
        **/
        void read_word_sym_table(const string& symFile);

        /** Clear out object except for symbol table. **/
        void clear();

        /** Returns whether there are no states. **/
        bool empty() const { return m_stateMap.empty(); }


        /** Get a pointer to the word symbol table. **/
        const SymbolTable& get_word_sym_table() const
            { return *m_symTable.get(); }

        /** Returns one above highest GMM index in graph. **/
        int get_gmm_count() const;

        /** Returns total number of states. **/
        int get_state_count() const { return m_stateMap.size(); }

        /** Returns index of start state. **/
        int get_start_state() const { return m_start; }

        /** Returns number of outgoing arcs for state @p stateIdx. **/
        int get_arc_count(int stateIdx) const
            {
            assert((stateIdx >= 0) && (stateIdx < (int) m_stateMap.size()));
            return get_max_arc_index(stateIdx) - get_min_arc_index(stateIdx);
            }

        /** Returns ref to @p arcIdx-th outgoing arc of state @p stateIdx. **/
        const Arc& get_arc(int stateIdx, int arcIdx) const
            {
            assert((stateIdx >= 0) && (stateIdx < (int) m_stateMap.size()));
            int minArcIdx = get_min_arc_index(stateIdx);
            assert((arcIdx >= 0) &&
                (arcIdx < get_max_arc_index(stateIdx) - minArcIdx));
            return m_arcList[minArcIdx + arcIdx];
            }

        /** Returns global index of an arc. **/
        int get_arc_index(const Arc& arc) const
            {
            assert((&arc >= &*m_arcList.begin()) &&
                (&arc < &*m_arcList.end()));
            return &arc - &*m_arcList.begin();
            }

        /** Returns ref to an arc given its global index. **/
        const Arc& get_arc(int arcIdx) const
            {
            assert((arcIdx >= 0) && (arcIdx < (int) m_arcList.size()));
            return m_arcList[arcIdx];
            }

        /** Returns src state of an arc. **/
        int get_src_state(const Arc& arc) const;

        /** Returns wheter state @p stateIdx is final state. **/
        bool is_final_state(int stateIdx) const
            { return m_finalCosts.find(stateIdx) != m_finalCosts.end(); }

        /** Returns final cost of state @p stateIdx or -g_zeroLogProb
        *   if not final.
        **/
        double get_final_cost(int stateIdx) const
            {
            map<int, float>::const_iterator lookup =
                m_finalCosts.find(stateIdx);
            return (lookup != m_finalCosts.end()) ? lookup->second :
                -g_zeroLogProb;
            }

        /** Returns number of final states; places list in @p stateList. **/
        int get_final_state_list(vector<int>& stateList) const
            {
            stateList.clear();
            for (map<int, float>::const_iterator elemPtr =
                m_finalCosts.begin(); elemPtr != m_finalCosts.end();
                ++elemPtr)
                stateList.push_back(elemPtr->first);
            sort(stateList.begin(), stateList.end());
            return stateList.size();
            }

    private:
        /** Returns index of first arc for state @p stateIdx. **/
        int get_min_arc_index(int stateIdx) const
            { return m_stateMap[stateIdx]; }

        /** Returns one past index of last arc for state @p stateIdx. **/
        int get_max_arc_index(int stateIdx) const
            { return (stateIdx + 1 < (int) m_stateMap.size()) ?
                m_stateMap[stateIdx + 1] : m_arcList.size(); }

    private:
        /** Map from words to integer indices. **/
        shared_ptr<SymbolTable> m_symTable;

        /** Index of start state. **/
        int m_start;

        /** Map from final states to their final costs. **/
        map<int, float> m_finalCosts;

        /** For each state, index of first arc in m_arcList.
        *   Assumes that arcs for each state are contiguous.
        **/
        vector<int> m_stateMap;

        /** List of arcs in graph, grouped by source state. **/
        vector<Arc> m_arcList;
    };


/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
*
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#endif


