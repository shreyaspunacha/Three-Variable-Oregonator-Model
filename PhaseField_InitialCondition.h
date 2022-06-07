double PhaseField_InitialCondition(int nx, int ny, int Circle_x, int Circle_y, int ObstacleRadius);
double PhaseFieldTimeEvolution(int nx, int ny, int nt, double dx, double dy, double dt, double zi,int Circle_x,int Circle_y,int ObstacleRadius);
double OregonatorTimeEvolution(int nx, int ny, double dx, double dy, double dt, int nt, double f, double q, double epsilon, double epsilonDash, int Circle_x, int Circle_y, int ObstacleRadius, double E, int PulseInterval, int savingInterval, int startTime, double Mu, double Mv, double Mw, int PulseStartTime);

