#pragma once

#include <iostream>
#include <map>
#include "Teuchos_RCP.hpp"
#ifdef HAVE_MPI
#include "Epetra_MpiComm.h"
#include <mpi.h>
#else
#include "Epetra_SerialComm.h"
#endif
#include "Teuchos_ParameterList.hpp"
#include "AnasaziTypes.hpp"
#include "AnasaziBasicEigenproblem.hpp"
#include "AnasaziGeneralizedDavidsonSolMgr.hpp"
#include "AnasaziEpetraAdapter.hpp"
#include "BelosSolverFactory.hpp"
#include "BelosSolverManager.hpp"
#include "BelosEpetraAdapter.hpp"
#include "Epetra_MultiVector.h"
#include "Epetra_Operator.h"
#include <cassert>

using std::map;

class AnasaziCnt{
public:
    /*
    Notes:
      solver manager needs a problem
      problem needs matrices
      the two above statements mean we can't constrcut them until solve time (consistent with Steven)
    */

    //We're solving   LHS*x = k*RHS*x
    Teuchos::RCP<Epetra_CrsMatrix> LHS;
    Teuchos::RCP<Epetra_CrsMatrix> RHS;
    Teuchos::RCP<Epetra_Operator>  pc;
    bool haspc=false;
    double keff;
    int niters;
    Teuchos::RCP<Epetra_Vector> x;
    Teuchos::ParameterList anasazi_db;
    //maybe some other things about a specific solver
};

class AnasaziStore {
public:
    AnasaziStore():
        cid(0)
    {}

    int new_data() {
        //Teuchos::ParameterList db
        things_[cid]=AnasaziCnt();
        //setup parameterlist with defaults
        //things_[cid].anasazi_db = Teuchos::sublist(db, "Anasazi");
        things_[cid].anasazi_db.set("Which", std::string("LM"));
        things_[cid].anasazi_db.get("Convergence Tolerance",1e-6);
        things_[cid].anasazi_db.get("Maximum Subspace Dimension",25);
        things_[cid].anasazi_db.get("Restart Dimension",5);
        things_[cid].anasazi_db.get("Maximum Restarts",100);
        things_[cid].anasazi_db.get("Initial Guess",std::string("User"));
        things_[cid].anasazi_db.get("Verbosity",Anasazi::Errors + Anasazi::Warnings);
        // + Anasazi::FinalSummary + Anasazi::TimingDetails

        cid++;
        return cid-1;
    }

    int delete_data(const int id){
        things_.erase(id);
        return 0;
    }

    int setMat_data(const int id, Teuchos::RCP<Epetra_CrsMatrix> LHS, Teuchos::RCP<Epetra_CrsMatrix> RHS) {
        things_[id].LHS=LHS;
        things_[id].RHS=RHS;
        return 0;
    }

    int setConvCrit_data(const int id, const double tol, const int maxit) {
        things_[id].anasazi_db.set("Convergence Tolerance", tol);
        things_[id].anasazi_db.set("Maximum Restarts", maxit);
        return 0;
    }

    int setX0_data(const int id, Teuchos::RCP<Epetra_Vector> x0) {
        things_[id].x=x0;
        return 0;
    }

    int setPC_data(const int id, Teuchos::RCP<Epetra_Operator> pc) {
        //if(things_[id].LHS->Comm().MyPID()==0) std::cout << pc->Label() << std::endl;
        things_[id].pc=pc;
        things_[id].haspc=true;
        return 0;
    }

    int solve(const int id) {
        Teuchos::RCP<Anasazi::BasicEigenproblem<double,Epetra_MultiVector,Epetra_Operator>> problem(
            new Anasazi::BasicEigenproblem<double,Epetra_MultiVector,Epetra_Operator>());
        problem->setA(things_[id].LHS);
        problem->setM(things_[id].RHS);
        if(things_[id].haspc) problem->setPrec(things_[id].pc);
        problem->setInitVec(things_[id].x);
        problem->setNEV(1);
        bool problem_set = problem->setProblem();
        assert(problem_set);

        Anasazi::GeneralizedDavidsonSolMgr<double,Epetra_MultiVector,Epetra_Operator> solver(
            problem, things_[id].anasazi_db);

        Anasazi::ReturnType returnval = solver.solve();
        assert(returnval==0);

        things_[id].niters=solver.getNumIters();

        // Extract solution
        Anasazi::Eigensolution<double,Epetra_MultiVector> solution =
            solver.getProblem().getSolution();
        Anasazi::Value<double> eval = (solution.Evals)[0];
        things_[id].keff = eval.realpart;
        things_[id].x->Update(1.0,*(solution.Evecs),0.0);
        //things_[id].x=Teuchos::rcp((*solution.Evecs)(0));

        return 0;
    }

    int getEigenvalue_data(const int id,double &keff) {
        keff=things_[id].keff;
        return 0;
    }

    int getIterations_data(const int id,int &niters) {
        niters=things_[id].niters;
        return 0;
    }

    int getResidual(const int id,double &resid) {
        Teuchos::RCP<Epetra_Vector> ltmp(new Epetra_Vector(*things_[id].x));
        Teuchos::RCP<Epetra_Vector> rtmp(new Epetra_Vector(*things_[id].x));
        things_[id].LHS->Multiply(false,*(things_[id].x),*ltmp);
        things_[id].RHS->Multiply(false,*(things_[id].x),*rtmp);
        double denom[1];
        (things_[id].x)->Norm2(denom);
        rtmp->Update(-1.,*(ltmp),things_[id].keff);
        double resids[1];
        rtmp->Norm2(resids);
        resid=resids[0]/denom[0];
        return 0;
    }

private:
        int cid;
        map<int, AnasaziCnt> things_;
};

class BelosCnt{
public:
    /*
    Notes:
      solver manager needs a problem
      problem needs matrices
      the two above statements mean we can't constrcut them until solve time (consistent with Steven)
    */
    Teuchos::RCP<Epetra_CrsMatrix> A;
    Teuchos::RCP<Epetra_Operator>  pc;
    bool haspc=false;
    int num_iters;
    Teuchos::RCP<Epetra_Vector> x;
    Teuchos::RCP<Epetra_Vector> b;
    Teuchos::ParameterList belos_db;
    //maybe some other things about a specific solver


    ~BelosCnt(){
        A  = Teuchos::null;
        pc = Teuchos::null;
        x  = Teuchos::null;
        b  = Teuchos::null;
    }
};

class BelosStore {
public:
    BelosStore():
        cid(0)
    {}

    int new_data() {
        //Teuchos::ParameterList db
        things_[cid]=BelosCnt();
        //setup parameterlist with defaults

        //things_[cid].anasazi_db = Teuchos::sublist(db, "Anasazi");
        //things_[cid].belos_db.get("belos_type","Pseudo Block GMRES");
        things_[cid].belos_db.get("Convergence Tolerance",1e-6);
        things_[cid].belos_db.get("Maximum Iterations",250);
        things_[cid].belos_db.get("Verbosity",0);  //Belos::Warnings+Belos::FinalSummary+Belos::StatusTestDetails
        things_[cid].belos_db.get("Output Frequency",1);
        things_[cid].belos_db.get("Implicit Residual Scaling","Norm of RHS");
        things_[cid].belos_db.get("Explicit Residual Scaling","Norm of RHS");

        cid++;
        return cid-1;
    }

    int delete_data(const int id) {
        things_.erase(id);
        return 0;
    }

    int setMat_data(const int id, Teuchos::RCP<Epetra_CrsMatrix> A) {
        things_[id].A=A;
        return 0;
    }

    int setConvCrit_data(const int id, const double tol, const int maxit) {
        things_[id].belos_db.set("Convergence Tolerance", tol);
        things_[id].belos_db.set("Maximum Iterations", maxit);
        return 0;
    }

    int setX0_data(const int id, Teuchos::RCP<Epetra_Vector> x0) {
        things_[id].x=x0;
        return 0;
    }

    int setb_data(const int id, Teuchos::RCP<Epetra_Vector> b) {
        things_[id].b=b;
        return 0;
    }

    int setPC_data(const int id, Teuchos::RCP<Epetra_Operator> pc) {
        things_[id].pc=pc;
        things_[id].haspc=true;
        return 0;
    }

    int solve(const int id) {
        Belos::SolverFactory<double,Epetra_MultiVector,Epetra_Operator> factory;

        std::string type ="Pseudo Block GMRES";
        Teuchos::RCP<Belos::SolverManager<double,Epetra_MultiVector,Epetra_Operator> > solver =
            factory.create(type,Teuchos::rcpFromRef(things_[id].belos_db));

        // Create linear problem
        Teuchos::RCP<Belos::LinearProblem<double,Epetra_MultiVector,Epetra_Operator> > problem(
                new Belos::LinearProblem<double,Epetra_MultiVector,Epetra_Operator>() );
        problem->setOperator(things_[id].A);
        problem->setLHS(things_[id].x);
        problem->setRHS(things_[id].b);
        if(things_[id].haspc) problem->setRightPrec(things_[id].pc);
        problem->setProblem();

        solver->setParameters(Teuchos::rcpFromRef(things_[id].belos_db));
        solver->setProblem(problem);

        Belos::ReturnType result = solver->solve();
        things_[id].num_iters = solver->getNumIters();

        return 0;
    }

    int getIterations_data(const int id,int &niter) {
        niter=things_[id].num_iters;
        return 0;
    }

    int getResidual(const int id,double &resid) {
        Teuchos::RCP<Epetra_Vector> rtmp(new Epetra_Vector(*things_[id].x));
        things_[id].A->Multiply(false,*(things_[id].x),*rtmp);
        rtmp->Update(-1.,*(things_[id].b),1.);
        double resids[1];
        rtmp->Norm2(resids);
        resid=resids[0];
        return 0;
    }

private:
        int cid;
        map<int, BelosCnt> things_;
};
