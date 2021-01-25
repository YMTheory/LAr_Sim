#ifndef LArGroupVelocity_h
#define LArGroupVelocity_h

class LArGroupVelocity
{
    public:
        LArGroupVelocity();
        ~LArGroupVelocity();

        static void setdata(double data) {m_data = data;}
        static void setdataerr(double dataerr) {m_data_err = dataerr;}
        static void setp0(double p0) {m_p0 = p0;}
        static double getp0()        {return m_p0;}
        static void setp1(double p1) {m_p1 = p1;}
        static double getp1()        {return m_p1;}
        static void setp2(double p2) {m_p2 = p2;}
        static double getp2()        {return m_p2;}

        static void Initialize();
        static void Calculate();

        //static void SetParameters();
        static double GetChi2();

    private:
        static double m_data;
        static double m_data_err;
        static double m_calc;
    
        static double m_p0;
        static double m_p1;
        static double m_p2;
        

};

#endif
