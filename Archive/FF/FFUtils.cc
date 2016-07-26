double Integrate(double (*F)(double), double minX, double maxX){
  return F(minX)*maxX;
}
