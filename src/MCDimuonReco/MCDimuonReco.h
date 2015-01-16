
class TreeHandler;
class MCDimuonReco : public TObject
{
  public:
  MCDimuonReco(TreeHandler& in);
  bool findHighPtDimuon();
  float getDimuonMass() const { return fmass; }
  private:
  TreeHandler& fT;
  float fmass;
  ClassDef(MCDimuonReco,1)
};
