
//  $Id: lab1.C,v 1.6 2009/09/17 14:47:08 stanchen Exp $


#include "util.H"
#include "front_end.H"


/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
*
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

void main_loop(const char** argv)
    {
    //  Process command line arguments.
    map<string, string> params;
    process_cmd_line(argv, params);

    //  Initialize front end.
    FrontEnd frontEnd(params);

    //  Main loop.
    ifstream audioStrm(
        get_required_string_param(params, "audio_file").c_str());
    ofstream featStrm(
        get_required_string_param(params, "feat_file").c_str());
    matrix<double> inAudio, feats;
    while (audioStrm.peek() != EOF)
        {
        string idStr = read_float_matrix(audioStrm, inAudio);
        cout << "Processing utterance ID: " << idStr << endl;
        frontEnd.get_feats(inAudio, feats);
        write_float_matrix(featStrm, feats, idStr);
        }
    audioStrm.close();
    featStrm.close();
    }

int main(int argc, const char** argv)
    {
    //  Catch all exceptions; print error message if exception caught.
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


