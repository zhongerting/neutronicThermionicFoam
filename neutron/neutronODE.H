#include "argList.H"
#include "IOmanip.H"
#include "ODESystem.H"
#include "ODESolver.H"

namespace Foam
{
    class neutronODE 
: 
    public ODESystem
{
    public:
    // Foam::scalar FissPower = 1.0;

    // 中子计算参数，包括缓发中子份额，衰变常数
    public:
    std::vector<Foam::scalar> aveTemperature = {0.0, 0.0, 0.0, 0.0};

    std::vector<Foam::scalar> betaFra = {2.618e-4, 1.737e-3, 1.555e-3, 3.133e-3, 9.123e-4, 3.330e-4};
    std::vector<Foam::scalar> betaCon = {1.271e-2, 3.174e-2, 1.160e-1, 3.110e-2, 1.397e+0, 3.872e+0};
    Foam::scalar betaTotal = 0.0079321;

    std::vector<Foam::scalar> decayFra = {0.01728, 0.01365, 0.01024, 0.02304};
    std::vector<Foam::scalar> decayCon = {1.00e-5, 3.00e-3, 1.00e-2, 4.00e-2};
    Foam::scalar decayTotal = 0.06421;

    Foam::scalar neulife = 0.2e-4;

    Foam::scalar react = 0.0;
    
    public:
    neutronODE()
    {
    }

    Foam::label nEqns() const
    {
        return 11;
    }

    // 中子方程，共11个方程：
    // 1(0)-->裂变中子方程
    // 2~7(1~6)-->缓发中子方程
    // 8~11(7~10)-->衰变热方程
    void derivatives
    (
        const scalar x, 
        const scalarField& y, 
        const label li,
        scalarField& dydt
    ) const
    {
        // 变量还原
        Foam::scalar fissPower = y[0];
        std::vector<Foam::scalar> betaI = {y[1], y[2], y[3], y[4], y[5], y[6]};
        std::vector<Foam::scalar> decayI = {y[7], y[8], y[9], y[10]};

        // 构建中间变量
        Foam::scalar fissPowerDt = y[0];
        std::vector<Foam::scalar> betaIDt = {y[1], y[2], y[3], y[4], y[5], y[6]};
        std::vector<Foam::scalar> decayIDt = {y[7], y[8], y[9], y[10]};

        // 构建方程
        Foam::scalar sumBetaCstPowerBetaI = 0.0;
        for (int i = 0; i < int(betaCon.size()); ++i)
        {
            sumBetaCstPowerBetaI += betaCon[i] * betaI[i];
        }

        // 中子裂变方程
        fissPowerDt = (react - betaTotal) * fissPower / neulife + sumBetaCstPowerBetaI;
        // 缓发中子方程
        for (int i = 0; i < int(betaCon.size()); ++i)
        {
            betaIDt[i] = betaFra[i] * fissPower / neulife - betaCon[i] * betaI[i];
        }
        // 衰变热方程
        for (int i = 0; i < int(decayCon.size()); ++i)
        {
            decayIDt[i] = decayFra[i] * fissPower - decayCon[i] * decayI[i];
        }
        
        
        // 变量传递
        dydt[0] = fissPowerDt;
        for (int i = 0; i < int(betaCon.size()); ++i)
        {
            dydt[i + 1] = betaIDt[i];
        }
        for (int i = 0; i < int(decayCon.size()); ++i)
        {
            dydt[i + 7] = decayIDt[i];
        }
    }

    void jacobian
    (
        const scalar x,
        const scalarField& y,
        const label li,
        scalarField& dfdx,
        scalarSquareMatrix& dfdy
    ) const
    {
        for (int i = 0; i < nEqns(); ++i)
        {
            dfdx[i] = 0.0;
        }
        dfdy(0, 0) = (react - betaTotal) / neulife;
        dfdy(0, 1) = betaCon[0];
        dfdy(0, 2) = betaCon[1];
        dfdy(0, 3) = betaCon[2];
        dfdy(0, 4) = betaCon[3];
        dfdy(0, 5) = betaCon[4];
        dfdy(0, 6) = betaCon[5];
        dfdy(0, 7) = 0.0;
        dfdy(0, 8) = 0.0;
        dfdy(0, 9) = 0.0;
        dfdy(0, 10) = 0.0;

        
        dfdy(1, 0) = betaFra[0] / neulife;
        dfdy(1, 1) = -betaCon[0];
        dfdy(1, 2) = 0.0;
        dfdy(1, 3) = 0.0;
        dfdy(1, 4) = 0.0;
        dfdy(1, 5) = 0.0;
        dfdy(1, 6) = 0.0;
        dfdy(1, 7) = 0.0;
        dfdy(1, 8) = 0.0;
        dfdy(1, 9) = 0.0;
        dfdy(1, 10) = 0.0;

        dfdy(2, 0) = betaFra[0] / neulife;
        dfdy(2, 1) = 0.0;
        dfdy(2, 2) = -betaCon[1];
        dfdy(2, 3) = 0.0;
        dfdy(2, 4) = 0.0;
        dfdy(2, 5) = 0.0;
        dfdy(2, 6) = 0.0;
        dfdy(2, 7) = 0.0;
        dfdy(2, 8) = 0.0;
        dfdy(2, 9) = 0.0;
        dfdy(2, 10) = 0.0;

        dfdy(3, 0) = betaFra[0] / neulife;
        dfdy(3, 1) = 0.0;
        dfdy(3, 2) = 0.0;
        dfdy(3, 3) = -betaCon[2];
        dfdy(3, 4) = 0.0;
        dfdy(3, 5) = 0.0;
        dfdy(3, 6) = 0.0;
        dfdy(3, 7) = 0.0;
        dfdy(3, 8) = 0.0;
        dfdy(3, 9) = 0.0;
        dfdy(3, 10) = 0.0;

        dfdy(4, 0) = betaFra[0] / neulife;
        dfdy(4, 1) = 0.0;
        dfdy(4, 2) = 0.0;
        dfdy(4, 3) = 0.0;
        dfdy(4, 4) = -betaCon[3];
        dfdy(4, 5) = 0.0;
        dfdy(4, 6) = 0.0;
        dfdy(4, 7) = 0.0;
        dfdy(4, 8) = 0.0;
        dfdy(4, 9) = 0.0;
        dfdy(4, 10) = 0.0;

        dfdy(5, 0) = betaFra[0] / neulife;
        dfdy(5, 1) = 0.0;
        dfdy(5, 2) = 0.0;
        dfdy(5, 3) = 0.0;
        dfdy(5, 4) = 0.0;
        dfdy(5, 5) = -betaCon[4];
        dfdy(5, 6) = 0.0;
        dfdy(5, 7) = 0.0;
        dfdy(5, 8) = 0.0;
        dfdy(5, 9) = 0.0;
        dfdy(5, 10) = 0.0;

        dfdy(6, 0) = betaFra[0] / neulife;
        dfdy(6, 1) = 0.0;
        dfdy(6, 2) = 0.0;
        dfdy(6, 3) = 0.0;
        dfdy(6, 4) = 0.0;
        dfdy(6, 5) = 0.0;
        dfdy(6, 6) = -betaCon[5];
        dfdy(6, 7) = 0.0;
        dfdy(6, 8) = 0.0;
        dfdy(6, 9) = 0.0;
        dfdy(6, 10) = 0.0;

        dfdy(7, 0) = decayFra[0];
        dfdy(7, 1) = 0.0;
        dfdy(7, 2) = 0.0;
        dfdy(7, 3) = 0.0;
        dfdy(7, 4) = 0.0;
        dfdy(7, 5) = 0.0;
        dfdy(7, 6) = 0.0;
        dfdy(7, 7) = -decayCon[0];
        dfdy(7, 8) = 0.0;
        dfdy(7, 9) = 0.0;
        dfdy(7, 10) = 0.0;

        dfdy(8, 0) = decayFra[0];
        dfdy(8, 1) = 0.0;
        dfdy(8, 2) = 0.0;
        dfdy(8, 3) = 0.0;
        dfdy(8, 4) = 0.0;
        dfdy(8, 5) = 0.0;
        dfdy(8, 6) = 0.0;
        dfdy(8, 7) = 0.0;
        dfdy(8, 8) = -decayCon[1];
        dfdy(8, 9) = 0.0;
        dfdy(8, 10) = 0.0;

        dfdy(9, 0) = decayFra[0];
        dfdy(9, 1) = 0.0;
        dfdy(9, 2) = 0.0;
        dfdy(9, 3) = 0.0;
        dfdy(9, 4) = 0.0;
        dfdy(9, 5) = 0.0;
        dfdy(9, 6) = 0.0;
        dfdy(9, 7) = 0.0;
        dfdy(9, 8) = 0.0;
        dfdy(9, 9) = -decayCon[2];
        dfdy(9, 10) = 0.0;

        dfdy(10, 0) = decayFra[0];
        dfdy(10, 1) = 0.0;
        dfdy(10, 2) = 0.0;
        dfdy(10, 3) = 0.0;
        dfdy(10, 4) = 0.0;
        dfdy(10, 5) = 0.0;
        dfdy(10, 6) = 0.0;
        dfdy(10, 7) = 0.0;
        dfdy(10, 8) = 0.0;
        dfdy(10, 9) = 0.0;
        dfdy(10, 10) = -decayCon[3];
    }
};
}