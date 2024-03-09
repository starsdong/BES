void fitCosTheta()
{
  TF1 *funCosTheta = new TF1("funCosTheta","[1]*((1-[0])+(3*[0]-1)*x*x)",-1,1);

  const *Char_t Name[3] = {"hPtCosTheta","hPtCosThetaRc","hPtCosThetaRcRes"};
  funCosTheta->SetParameters(rho3/3.);
}
