void a_half(Double_t BL = 0.13, Double_t r_half = 2.3)  // BL in MeV, r_half in fm
{
  const Double_t hbarc = 197.3;  // MeV*fm
  const Double_t m_d = 1875.6;  // MeV - deuteron mass
  const Double_t m_l = 1115.7;  // MeV - Lambda mass
  Double_t mu = m_d*m_l/(m_d + m_l);  // reduced mass - MeV

  Double_t gamma = TMath::Sqrt(2*mu*BL);  // MeV

  double inv_a_half = gamma - 0.5 * r_half * gamma * gamma / hbarc;  // MeV
  //  double inv_a_half = gamma;

  double a_half = 1/inv_a_half * hbarc;

  std::cout << " gamma = " << gamma << "\t a_half = " << a_half << std::endl;
  return a_half;

}

void BL(Double_t a_half = 10.0, Double_t r_half = 2.3)  // BL in MeV, r_half in fm, a_half in fm
{
  const Double_t hbarc = 197.3;  // MeV*fm
  const Double_t m_d = 1875.6;  // MeV - deuteron mass
  const Double_t m_l = 1115.7;  // MeV - Lambda mass
  Double_t mu = m_d*m_l/(m_d + m_l);  // reduced mass - MeV


  Double_t gamma = (1 - sqrt(1-2*r_half/a_half))/r_half *hbarc;  // MeV

  double BL = gamma * gamma / 2./ mu;

  std::cout << " gamma = " << gamma << "\t BL = " << BL << std::endl;
  return BL;

}
