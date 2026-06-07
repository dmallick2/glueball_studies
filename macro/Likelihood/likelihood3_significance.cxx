#include <iostream>
#include <cmath>

#include "TH1D.h"
#include "TRandom3.h"
#include "TMath.h"

//// In this code we calculating the significance of the signal by comparing the likelihood of the signal+background model with the background-only model. We are doing this by scanning the parameter space for both cases and finding the best likelihood for each case. Then we can use these likelihood values to compute the significance using the formula: significance = sqrt(2 * (logL_s+b - logL_b)).
///// *** Can change the counts of signal and background to see how the significance changes.***

//// In this code we will also calculate the p-value of the fit using the likelihood ratio test. p value is the probability of obtaining a test statistic (q) equal to or more than the observed value under the null hypothesis (background only). Suppose p value is large, then it means that the observed data is consistent with the background only hypothesis and we do not have evidence for the signal. If p-value is small, then means that it is unlikely to observe a signal event from background fluctuation and we have evidence for the signal. In particle physics, we usually require a p-value of less than 2.87e-7 (corresponding to a significance of 5 sigma) to claim a discovery.

using namespace std;

//--------------------------------------------------
// Signal shape
//--------------------------------------------------

double gaussian(double x, double mean, double sigma)
{
    return TMath::Gaus(x, mean, sigma, true);
}

//--------------------------------------------------
// Background shape
//--------------------------------------------------

double background(double x)
{
    return (0.1 + 0.02 * x) / 2.0; // Factor 2.0 is for normalization over the range 0-10
}

//--------------------------------------------------
// Total model
//--------------------------------------------------

double totalFunction(double *x, double *par)
{
    double xx = x[0];

    double Ns = par[0];
    double Nb = par[1];

    double mean = 5.0;
    double sigma = 0.3;

    double s = gaussian(xx, mean, sigma);
    double b = background(xx);

    // expected counts per bin
    double bw = 0.1;

    return (Ns * s + Nb * b) * bw;
}

//--------------------------------------------------
// Main function
//--------------------------------------------------

void likelihood3_significance()
{
    //--------------------------------------------------
    // Create toy histogram
    //--------------------------------------------------

    TH1D *h = new TH1D("h", "Toy data", 100, 0, 10);

    TRandom3 r(0);

    //--------------------------------------------------
    // Generate signal
    //--------------------------------------------------

    for (int i = 0; i < 1000; i++)
    {
        double x = r.Gaus(5.0, 0.3);
        h->Fill(x);
    }

    //--------------------------------------------------
    // Generate background
    //--------------------------------------------------

    TF1 *fbkg = new TF1("fbkg", "0.1 + 0.02*x", 0, 10);

    for (int i = 0; i < 9000; i++)
    {
        double x = fbkg->GetRandom();
        h->Fill(x);
    }

    //--------------------------------------------------
    // Fixed shape parameters
    //--------------------------------------------------

    double mean = 5.0;
    double sigma = 0.3;

    //--------------------------------------------------
    // Variables to store best fit
    //--------------------------------------------------

    double bestLogL = -1e30;

    double bestNs = 0;
    double bestNb = 0;

    //--------------------------------------------------
    // Scan parameter space
    //--------------------------------------------------

    for (double Ns = 950; Ns <= 1030; Ns += 1)
    {
        for (double Nb = 8850; Nb <= 10000; Nb += 1)
        {
            double logL = 0;

            //--------------------------------------------------
            // Loop over histogram bins
            //--------------------------------------------------

            for (int i = 1; i <= h->GetNbinsX(); i++)
            {
                double x = h->GetBinCenter(i);
                double bw = h->GetBinWidth(i);

                // observed counts
                double n = h->GetBinContent(i);

                //--------------------------------------------------
                // Expected counts
                //--------------------------------------------------

                double s = gaussian(x, mean, sigma);

                double b = background(x);

                double mu = (Ns * s + Nb * b) * bw;

                //--------------------------------------------------
                // Poisson log likelihood
                //--------------------------------------------------

                if (mu > 0)
                {
                    logL += n * log(mu) - mu;
                }
            }

            //--------------------------------------------------
            // Check if likelihood improved
            //--------------------------------------------------

            if (logL > bestLogL)
            {
                bestLogL = logL;

                bestNs = Ns;
                bestNb = Nb;
            }
        }
    }

    //--------------------------------------------------
    // Print best fit
    //--------------------------------------------------

    cout << "Best fit results" << endl;

    cout << "Signal yield     = " << bestNs << endl;

    cout << "Background yield = " << bestNb << endl;

    cout << "Maximum logL     = " << bestLogL << endl;

    //--------------------------------------------------
    // Draw histogram and the best fit
    //--------------------------------------------------
    h->SetLineColor(kBlack);
    h->Draw();
    TF1 *f = new TF1("f", totalFunction, 0, 10, 2);
    f->SetParameter(0, bestNs);
    f->SetParameter(1, bestNb);
    f->SetLineStyle(2);
    f->Draw("same");

    //--------------------------------------------------
    // To compute the significance we compare the signal+background model with background only model.
    //--------------------------------------------------

    //--------------------------------------------------
    // BACKGROUND ONLY FIT
    //--------------------------------------------------

    double bestLogL_B = -1e30;

    double bestNb_B = 0;

    for (double Nb = 8850; Nb <= 10000; Nb += 1)
    {
        double logL = 0;

        for (int i = 1; i <= h->GetNbinsX(); i++)
        {
            double x = h->GetBinCenter(i);
            double bw = h->GetBinWidth(i);

            double n = h->GetBinContent(i);

            //--------------------------------------
            // signal fixed to ZERO
            //--------------------------------------

            double b = background(x);

            double mu = (Nb * b) * bw;

            if (mu > 0)
            {
                logL += n * log(mu) - mu;
            }
        }

        //--------------------------------------
        // keep best background-only likelihood
        //--------------------------------------

        if (logL > bestLogL_B)
        {
            bestLogL_B = logL;
            bestNb_B = Nb;
        }
    }

    //--------------------------------------------------
    // LIKELIHOOD RATIO SIGNIFICANCE
    //--------------------------------------------------

    double q = 2.0 * (bestLogL - bestLogL_B);

    double significance = sqrt(q);

    //--------------------------------------------------
    // p-value from Gaussian tail
    //--------------------------------------------------

    double pvalue = 0.5 * TMath::Erfc(significance / sqrt(2.0)); // Two-sided p-value for Gaussian distribution

    cout << endl;
    cout << "===== SIGNIFICANCE =====" << endl;
    cout << "Best logL (S+B) = " << bestLogL << endl;
    cout << "Best logL (B)   = " << bestLogL_B << endl;
    cout << "q = " << q << endl;
    cout << "Significance = " << significance << " sigma" << endl;
    cout << "p-value = " << pvalue << endl;

    //--------------------------------------------------
    // Same approach with ROOT's option "L" for likelihood fit
    //--------------------------------------------------

    TF1 *fROOT = new TF1("fROOT", totalFunction, 0, 10, 2);
    fROOT->SetParameter(0, 1000);
    fROOT->SetParameter(1, 9000);
    TFitResultPtr fitResult = h->Fit(fROOT, "LS0");
    fitResult->Print("V"); // This will print the fit result and the correlation matrix. The "V" option is for verbose output. Correlation matrix tells us about how much the parameters are correlated with each other. The off diagonal elements range from -1 to 1, where -1 means perfect anti-correlation, 0 means no correlation, and 1 means perfect correlation. If the parameters are highly correlated, it can make the fit less stable and the uncertainties on the parameters can be larger. Therefore if the value is around 0.5 or less, then it means that the parameters are not highly correlated and the fit is more stable.

    TF1 *fROOT_B = new TF1("fROOT_B", totalFunction, 0, 10, 2);
    fROOT_B->FixParameter(0, 0); // Signal fixed to zero for background-only fit
    fROOT_B->SetParameter(1, 9000);
    TFitResultPtr fitResult_B = h->Fit(fROOT_B, "LS0");

    double minLogL_SB = fitResult->MinFcnValue(); // ROOT returns -logL
    double minLogL_B = fitResult_B->MinFcnValue();
    double q_ROOT = 2.0 * (minLogL_B - minLogL_SB); // Note the order of logL_B and logL_SB for ROOT as it returns -logL
    double significance_ROOT = sqrt(q_ROOT);
    double pvalue_ROOT = 0.5 * TMath::Erfc(significance_ROOT / sqrt(2.0)); // Two-sided p-value for Gaussian distribution

    cout << endl;
    cout << "===== ROOT OPTION L FIT =====" << endl;

    cout << "Signal     = " << fROOT->GetParameter(0) << endl;
    cout << "Background = " << fROOT->GetParameter(1) << endl;
    cout << "Minimum -logL = " << fitResult->MinFcnValue() / 2.0 << endl;
    cout << "q (ROOT) = " << q_ROOT << endl;
    cout << "Significance (ROOT) = " << significance_ROOT << " sigma" << endl;
    cout << "p-value (ROOT) = " << pvalue_ROOT << endl;
}

////=============Output=============
// Correlation Matrix:

//             	          p0          p1
// p0          	           1    -0.25052
// p1          	    -0.25052           1
// The output shows that signal and background are slightly anti-correlated, which is expected because an increase in signal yield can be partially compensated by a decrease in background yield to fit the data. However, since the correlation is not very strong (around -0.25), it indicates that the fit is reasonably stable.
