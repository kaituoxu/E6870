
//  $Id: lab1_dtw.C,v 1.10 2009/09/18 02:12:13 stanchen Exp $


#include "util.H"



/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
*
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

double compute_distance(const matrix<double>& matHyp,
    const matrix<double>& matTempl)
    {
    double dist = 0.0;

    //  BEGIN_LAB
    //
    //  Input:
    //      "matHyp", matrix containing feature vectors for the current
    //      test utterance, and "matTempl", matrix containing feature
    //      vectors for the current template utterance being compared
    //      with.  These two matrices are guaranteed to have the same
    //      number of columns (but not the same number of frames).
    //
    //      For a matrix "mat", "mat.size1()" returns the number of
    //      frames, "mat.size2()" returns the number of dimensions in
    //      each feature vector, and "mat(frmIdx, dimIdx)" is the
    //      syntax for accessing elements.  Indices are numbered from 0.
    //
    //  Output:
    //      Set "dist" to the total distance between these two utterances.

    int inRow = matHyp.size1();
    int inCol = matHyp.size2();
    int outRow = matTempl.size1();
    int outCol = matTempl.size2();
    int i,j;
    matrix<double> distance;
    matrix<double> output;

    distance.resize(inRow, outRow);
    output.resize(inRow + 1, outRow + 1);

    for(i = 1;i <= inRow;i++){
      output(i, 0) = 1000000.0;
    }
    for(j = 1;j <= outRow;j++){
      output(0, j) = 1000000.0;
    }
    output(0, 0) = 0;

    for(i = 0;i < inRow;i++){
        for(j = 0;j < outRow;j++){
            double sum = 0;
            for(int idx = 0;idx < inCol;idx++){
	      sum = sum + pow((matHyp(i, idx) - matTempl(j,idx)),2);
            }
            distance(i,j) = sqrt(sum);
        }
    }

    for(i = 1;i <= inRow;i++)
        for(j = 1;j <= outRow;j++)
	  output(i,j) = min(min(output(i - 1,j - 1),output(i,j - 1)),output(i - 1,j)) + distance(i - 1,j - 1);

    dist = output(inRow, outRow);


    //  END_LAB

    return dist;
    }


/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
*
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

void main_loop(const char** argv)
    {
    //  Process command line arguments.
    map<string, string> params;
    process_cmd_line(argv, params);
    bool verbose = get_bool_param(params, "verbose");

    //  Load feature files for templates.
    //  Get template label from matrix name.
    vector<string> templateLabelList;
    vector<matrix<double> > templateMatList;
    ifstream templateStrm(
        get_required_string_param(params, "template_file").c_str());
    while (templateStrm.peek() != EOF)
        {
        templateMatList.push_back(matrix<double>());
        string labelStr = read_float_matrix(templateStrm,
            templateMatList.back());
        templateLabelList.push_back(labelStr);
        }
    templateStrm.close();
    if (templateMatList.empty())
        throw runtime_error("No templates supplied.");

    //  Load correct label for each feature file, if present.
    vector<string> featLabelList;
    string labelFile = get_string_param(params, "feat_label_list");
    if (!labelFile.empty())
        read_string_list(labelFile, featLabelList);

    //  The main loop.
    ifstream featStrm(
        get_required_string_param(params, "feat_file").c_str());
    matrix<double> feats;
    unsigned templCnt = templateLabelList.size();
    unsigned uttCnt = 0;
    unsigned correctCnt = 0;
    while (featStrm.peek() != EOF)
        {
        int uttIdx = uttCnt++;
        string idStr = read_float_matrix(featStrm, feats);

        //  Find closest template.
        int bestTempl = -1;
        double bestScore = DBL_MAX;
        for (unsigned templIdx = 0; templIdx < templCnt; ++templIdx)
            {
            if (feats.size2() != templateMatList[templIdx].size2())
                throw runtime_error("Mismatch in test/template feature dim.");
            double curScore = compute_distance(feats,
                templateMatList[templIdx]);
            if (verbose)
                cout << format("  %s: %.3f") % templateLabelList[templIdx] %
                    curScore << endl;
            if (curScore < bestScore)
                bestScore = curScore, bestTempl = templIdx;
            }
        if (bestTempl < 0)
            throw runtime_error("No alignments found in DTW.");

        string hypLabel = (bestTempl >= 0) ? templateLabelList[bestTempl] : "";
        if (!featLabelList.empty())
            {
            //  If have reference labels, print ref and hyp classes.
            if (uttIdx >= (int) featLabelList.size())
                throw runtime_error("Mismatch in number of utterances "
                    "and labels.");
            string refLabel = featLabelList[uttIdx];
            cout << format("Reference: %s, Hyp: %s, Correct: %d") %
                refLabel % hypLabel % (hypLabel == refLabel) << endl;
            if (hypLabel == refLabel)
                ++correctCnt;
            }
        else
            //  If don't have reference labels, just print hyp class.
            cout << hypLabel << " (" << idStr << ")" << endl;
        }
    featStrm.close();
    if (!featLabelList.empty())
        {
        //  If have reference labels, print accuracy.
        unsigned errCnt = uttCnt - correctCnt;
        cout << format("Accuracy: %.2f%% (%d/%d), Error rate: %.2f%% (%d/%d)")
            % (100.0 * correctCnt / uttCnt) % correctCnt % uttCnt %
            (100.0 * errCnt / uttCnt) % errCnt % uttCnt << endl;
        }
    }

int main(int argc, const char** argv)
    {
    try
        {
        main_loop(argv);
        }
    catch (exception& xc)
        {
        cerr << "Error: " << xc.what() << endl;
        return -1;
        }
    return 0;
    }


/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
*
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */


