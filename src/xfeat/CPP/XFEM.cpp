int EqNum(int * JJglob, int Njoint) {
	int Ndf = 0;
	for (int i = 1; i < Njoint + 1; i++)
		for (int k = 0; k < 3; k++)
			if (JJglob[(3 * i) + k] >= 0) {
				JJglob[(3 * i) + k] = Ndf;
				++Ndf;
			}
	return Ndf;
}

void JacInv() {

	// FORMATION OF THE JACOBIAN INVERSION MATRIX

	int i, j;
	double C[4][4];

	C[1][1] = Ajac[2][2] * Ajac[3][3] - Ajac[2][3] * Ajac[3][2];
	C[1][2] = -(Ajac[2][1] * Ajac[3][3] - Ajac[2][3] * Ajac[3][1]);
	C[1][3] = Ajac[2][1] * Ajac[3][2] - Ajac[2][2] * Ajac[3][1];
	C[2][1] = -(Ajac[1][2] * Ajac[3][3] - Ajac[1][3] * Ajac[3][2]);
	C[2][2] = Ajac[1][1] * Ajac[3][3] - Ajac[1][3] * Ajac[3][1];
	C[2][3] = -(Ajac[1][1] * Ajac[3][2] - Ajac[1][2] * Ajac[3][1]);
	C[3][1] = Ajac[1][2] * Ajac[2][3] - Ajac[1][3] * Ajac[2][2];
	C[3][2] = -(Ajac[1][1] * Ajac[2][3] - Ajac[1][3] * Ajac[2][1]);
	C[3][3] = Ajac[1][1] * Ajac[2][2] - Ajac[1][2] * Ajac[2][1];

	Det = 0.0;

	for (i = 1; i < 4; i++) {
		Det = Det + Ajac[1][i] * C[1][i];
	}

	for (i = 1; i < 4; i++) {
		for (j = 1; j < 4; j++) {
			AjacInv[i][j] = C[j][i] / Det;
		}
	}
}

void Assem_force(int Iel) {
	int i, icon, iicon, ii, iii;

	for (i = 0; i < 8; i++) {
		icon = LOTOGO[8*Iel + i];
		for (ii = 0; ii < 3; ii++) {
			iicon = JJglob[3*icon + ii];
			iii = 3*i + ii;
			if (iicon >= 0)
				Fglob[iicon] += ENFLOC[iii+1];
		}
	}
}

void AssemDislocationForce(int Iel) {
	int i, icon, iicon, ii, iii;

	for (i = 0; i < 8; i++) {
		icon = LOTOGO[8*Iel + i];
		for (ii = 0; ii < 3; ii++) {
			iicon = JJglob[(3 * icon) + ii];
			iii = 3*i + ii;
			if (iicon >= 0)
				Fglob[iicon] += DISLOCATIONF[iii+1];
		}
	}
}

void ENFORCEDFORCE(double XLOC[NODEN], double YLOC[NODEN], double ZLOC[NODEN],
		int IEL) {

	double A[NDS][10], B[10][10], C[10][NDLOC];
	double BB[NDS][NDLOC], BBTRAN[NDLOC][NDS];

	double SND[ND][NODEN];

	double BBTRANE[NDLOC][NDS], AB[NDS][10];

	double SHI, NEW, KSHI;
	double DIR1[NODEN], DIR2[NODEN], DIR3[NODEN], COR[NODEN][ND];

	DIR1[1] = DIR2[4] = DIR3[3] = -1.0;
	DIR1[2] = DIR2[8] = DIR3[4] = -1.0;
	DIR1[3] = DIR2[3] = DIR3[7] = 1.0;
	DIR1[4] = DIR2[7] = DIR3[8] = 1.0;
	DIR1[5] = DIR2[1] = DIR3[2] = -1.0;
	DIR1[6] = DIR2[5] = DIR3[1] = -1.0;
	DIR1[7] = DIR2[2] = DIR3[5] = 1.0;
	DIR1[8] = DIR2[6] = DIR3[6] = 1.0;

	for (int i = 0; i < NDS; i++)
		for (int j = 0; j < 10; j++)
			A[i][j] = 0.0;

	A[1][1] = A[2][5] = A[3][9] = 1.0;
	A[4][2] = A[4][4] = 1.0;
	A[5][6] = A[5][8] = 1.0;
	A[6][3] = A[6][7] = 1.0;

	for (int i = 0; i < NDLOC; i++)
		ENFLOC[i] = 0.0;

	for (int III = 1; III < 3; III++)        		// GAUSS INTEGRATION LOOP
		for (int JJJ = 1; JJJ < 3; JJJ++)
			for (int PPP = 1; PPP < 3; PPP++) {
				if (III == 1)
					SHI = -0.577350269189626;
				if (III == 2)
					SHI = 0.577350269189626;

				if (JJJ == 1)
					NEW = -0.577350269189626;
				if (JJJ == 2)
					NEW = 0.577350269189626;

				if (PPP == 1)
					KSHI = -0.577350269189626;
				if (PPP == 2)
					KSHI = 0.577350269189626;

				for (int i = 0; i < ND; i++)
					for (int j = 0; j < NODEN; j++)
						SND[i][j] = 0.0;

				for (int j = 1; j < NODEN; j++) {
					//Order changed to match Abaqus numbering
					SND[2][j] = (1.0 / 8.0) * DIR1[j] * (1.0 + NEW * DIR2[j])
							* (1.0 + KSHI * DIR3[j]);
					SND[1][j] = (1.0 / 8.0) * DIR2[j] * (1.0 + SHI * DIR1[j])
							* (1.0 + KSHI * DIR3[j]);
					SND[3][j] = (1.0 / 8.0) * DIR3[j] * (1.0 + NEW * DIR2[j])
							* (1.0 + SHI * DIR1[j]);
				}

				for (int i = 0; i < NODEN; i++)
					for (int j = 0; j < ND; j++)
						COR[i][j] = 0.0;

				for (int i = 1; i < NODEN; i++) {
					COR[i][1] = XLOC[i];
					COR[i][2] = YLOC[i];
					COR[i][3] = ZLOC[i];
				}

				for (int i = 0; i < ND; i++)
					for (int j = 0; j < ND; j++) {
						Ajac[i][j] = 0.0;
						AjacInv[i][j] = 0.0;
					}

				for (int i = 1; i < ND; i++)
					for (int j = 1; j < ND; j++)
						for (int k = 1; k < NODEN; k++)
							Ajac[i][j] = Ajac[i][j] + SND[i][k] * COR[k][j];

				JacInv();

				for (int i = 0; i < 10; i++)
					for (int j = 0; j < 10; j++)
						B[i][j] = 0.0;

				for (int j = 1; j < ND; j++)
					for (int l = 1; l < ND; l++) {
						B[j][l] = AjacInv[j][l];
						B[j + 3][l + 3] = AjacInv[j][l];
						B[j + 6][l + 6] = AjacInv[j][l];
					}

				for (int i = 0; i < 10; i++)
					for (int j = 0; j < NDLOC; j++)
						C[i][j] = 0.0;

				for (int i = 1; i < 4; i++)
					for (int j = 1; j < 24; j = j + 3)
						C[i][j] = SND[i][(j + 2) / 3];

				for (int i = 4; i < 7; i++)
					for (int j = 2; j < 24; j = j + 3)
						C[i][j] = SND[i - 3][(j + 1) / 3];

				for (int i = 7; i < 10; i++)
					for (int j = 3; j < 25; j = j + 3)
						C[i][j] = SND[i - 6][j / 3];

				for (int i = 0; i < NDS; i++)
					for (int j = 0; j < 10; j++)
						AB[i][j] = 0.0;

				for (int i = 1; i < NDS; i++)
					for (int j = 1; j < 10; j++)
						for (int k = 1; k < 10; k++)
							AB[i][j] = AB[i][j] + A[i][k] * B[k][j];

				for (int i = 0; i < NDS; i++)
					for (int j = 0; j < NDLOC; j++)
						BB[i][j] = 0.0;

				for (int i = 1; i < NDS; i++)
					for (int j = 1; j < NDLOC; j++)
						for (int k = 1; k < 10; k++)
							BB[i][j] = BB[i][j] + AB[i][k] * C[k][j];

				for (int i = 0; i < NDLOC; i++)
					for (int j = 0; j < NDS; j++) {
						BBTRAN[i][j] = 0.0;
						BBTRANE[i][j] = 0.0;
					}

				for (int i = 1; i < NDS; i++)
					for (int j = 1; j < NDLOC; j++)
						BBTRAN[j][i] = BB[i][j];

				for (int i = 1; i < NDLOC; i++)
					for (int j = 1; j < NDS; j++)
						for (int k = 1; k < NDS; k++)
							BBTRANE[i][j] = BBTRANE[i][j]
									+ BBTRAN[i][k] * EE3D[k][j];

				double ENDP[NDLOC], BBENDP[NDS];

				for (int i = 1; i < NODEN; i++) {
					ENDP[(3 * i - 2)] = ENFRDISPglob[3
							* LOTOGO[8 * IEL + (i - 1)] + (1 - 1)];
					ENDP[(3 * i - 1)] = ENFRDISPglob[3
							* LOTOGO[8 * IEL + (i - 1)] + (2 - 1)];
					ENDP[3 * i] = ENFRDISPglob[3
					        * LOTOGO[8 * IEL + (i - 1)] + (3 - 1)];
				}

				for (int i = 1; i < NDS; i++)
					BBENDP[i] = 0.0;

				for (int i = 1; i < NDS; i++)
					for (int k = 1; k < NDLOC; k++)
						BBENDP[i] = BBENDP[i] + BB[i][k] * ENDP[k];

				for (int i = 1; i < NDLOC; i++)
					for (int k = 1; k < NDS; k++)
						ENFLOC[i] = ENFLOC[i]
								- (BBTRANE[i][k] * BBENDP[k] * Det);
			}
}

void DISLOCATIONFORCE(double XLOC[NODEN], double YLOC[NODEN],
		double ZLOC[NODEN], int IEL, double bv, double et[ND]) {

	double A[NDS][10], B[10][10], C[10][NDLOC];
	double BB[NDS][NDLOC], BBTRAN[NDLOC][NDS];

	double SND[ND][NODEN];

	double BBTRANE[NDLOC][NDS], AB[NDS][10];

	double SHI, NEW, KSHI;
	double DIR1[NODEN], DIR2[NODEN], DIR3[NODEN], COR[NODEN][ND];

	DIR1[1] = DIR2[4] = DIR3[3] = -1.0;
	DIR1[2] = DIR2[8] = DIR3[4] = -1.0;
	DIR1[3] = DIR2[3] = DIR3[7] = 1.0;
	DIR1[4] = DIR2[7] = DIR3[8] = 1.0;
	DIR1[5] = DIR2[1] = DIR3[2] = -1.0;
	DIR1[6] = DIR2[5] = DIR3[1] = -1.0;
	DIR1[7] = DIR2[2] = DIR3[5] = 1.0;
	DIR1[8] = DIR2[6] = DIR3[6] = 1.0;

	for (int i = 0; i < NDS; i++)
		for (int j = 0; j < 10; j++)
			A[i][j] = 0.0;

	A[1][1] = A[2][5] = A[3][9] = 1.0;
	A[4][2] = A[4][4] = 1.0;
	A[5][6] = A[5][8] = 1.0;
	A[6][3] = A[6][7] = 1.0;

	for (int i = 0; i < NDLOC; i++)
		DISLOCATIONF[i] = 0.0;

	for (int III = 1; III < 3; III++)        		// GUASS INTEGRATION LOOP
		for (int JJJ = 1; JJJ < 3; JJJ++)
			for (int PPP = 1; PPP < 3; PPP++) {
				if (III == 1)
					SHI = -0.577350269189626;
				if (III == 2)
					SHI = 0.577350269189626;

				if (JJJ == 1)
					NEW = -0.577350269189626;
				if (JJJ == 2)
					NEW = 0.577350269189626;

				if (PPP == 1)
					KSHI = -0.577350269189626;
				if (PPP == 2)
					KSHI = 0.577350269189626;

				for (int i = 0; i < ND; i++)
					for (int j = 0; j < NODEN; j++)
						SND[i][j] = 0.0;

				for (int j = 1; j < NODEN; j++) {
					//Order changed to match Abaqus numbering
					SND[2][j] = (1.0 / 8.0) * DIR1[j] * (1.0 + NEW * DIR2[j])
							* (1.0 + KSHI * DIR3[j]);
					SND[1][j] = (1.0 / 8.0) * DIR2[j] * (1.0 + SHI * DIR1[j])
							* (1.0 + KSHI * DIR3[j]);
					SND[3][j] = (1.0 / 8.0) * DIR3[j] * (1.0 + NEW * DIR2[j])
							* (1.0 + SHI * DIR1[j]);
				}

				for (int i = 0; i < NODEN; i++)
					for (int j = 0; j < ND; j++)
						COR[i][j] = 0.0;

				for (int i = 1; i < NODEN; i++) {
					COR[i][1] = XLOC[i];
					COR[i][2] = YLOC[i];
					COR[i][3] = ZLOC[i];
				}

				for (int i = 0; i < ND; i++)
					for (int j = 0; j < ND; j++) {
						Ajac[i][j] = 0.0;
						AjacInv[i][j] = 0.0;
					}

				for (int i = 1; i < ND; i++)
					for (int j = 1; j < ND; j++)
						for (int k = 1; k < NODEN; k++)
							Ajac[i][j] = Ajac[i][j] + SND[i][k] * COR[k][j];

				JacInv();

				double GND[ND], H;

				for (int i = 0; i < ND; i++)
					GND[i] = 0.0;

				for (int node = 1; node < NODEN; node++) {
					for (int i = 1; i < ND; i++)
						for (int j = 1; j < ND; j++) {
							if (YLOC[node] < 0.0)
								H = -0.5;
							else
								H = 0.5;

							GND[i] = GND[i] + AjacInv[i][j] * SND[j][node] * H;
						}
				}

				for (int i = 0; i < 10; i++)
					for (int j = 0; j < 10; j++)
						B[i][j] = 0.0;

				for (int j = 1; j < ND; j++)
					for (int l = 1; l < ND; l++) {
						B[j][l] = AjacInv[j][l];
						B[j + 3][l + 3] = AjacInv[j][l];
						B[j + 6][l + 6] = AjacInv[j][l];
					}

				for (int i = 0; i < 10; i++)
					for (int j = 0; j < NDLOC; j++)
						C[i][j] = 0.0;

				for (int i = 1; i < 4; i++)
					for (int j = 1; j < 24; j = j + 3)
						C[i][j] = SND[i][(j + 2) / 3];

				for (int i = 4; i < 7; i++)
					for (int j = 2; j < 24; j = j + 3)
						C[i][j] = SND[i - 3][(j + 1) / 3];

				for (int i = 7; i < 10; i++)
					for (int j = 3; j < 25; j = j + 3)
						C[i][j] = SND[i - 6][j / 3];

				for (int i = 0; i < NDS; i++)
					for (int j = 0; j < 10; j++)
						AB[i][j] = 0.0;

				for (int i = 1; i < NDS; i++)
					for (int j = 1; j < 10; j++)
						for (int k = 1; k < 10; k++)
							AB[i][j] = AB[i][j] + A[i][k] * B[k][j];

				for (int i = 0; i < NDS; i++)
					for (int j = 0; j < NDLOC; j++)
						BB[i][j] = 0.0;

				for (int i = 1; i < NDS; i++)
					for (int j = 1; j < NDLOC; j++)
						for (int k = 1; k < 10; k++)
							BB[i][j] = BB[i][j] + AB[i][k] * C[k][j];

				for (int i = 0; i < NDLOC; i++)
					for (int j = 0; j < NDS; j++) {
						BBTRAN[i][j] = 0.0;
						BBTRANE[i][j] = 0.0;
					}

				for (int i = 1; i < NDS; i++)
					for (int j = 1; j < NDLOC; j++)
						BBTRAN[j][i] = BB[i][j];

				for (int i = 1; i < NDLOC; i++)
					for (int j = 1; j < NDS; j++)
						for (int k = 1; k < NDS; k++)
							BBTRANE[i][j] = BBTRANE[i][j]
									+ BBTRAN[i][k] * EE3D[k][j];

				double Dalpha[NDS], ex[ND], ey[ND], ez[ND];
				double ext, eyt, ezt;

				ext = eyt = ezt = 0.0;

				ex[1] = 1.0;
				ex[2] = 0.0;
				ex[3] = 0.0;

				ey[1] = 0.0;
				ey[2] = 1.0;
				ey[3] = 0.0;

				ez[1] = 0.0;
				ez[2] = 0.0;
				ez[3] = 1.0;

				for (int i = 1; i < ND; i++) {
					ext = ext + ex[i] * et[i];
					eyt = eyt + ey[i] * et[i];
					ezt = ezt + ez[i] * et[i];
				}

				Dalpha[1] = bv * ext * GND[1];
				Dalpha[2] = bv * eyt * GND[2];
				Dalpha[3] = bv * ezt * GND[3];
				Dalpha[4] = bv * (eyt * GND[1] + ext * GND[2]);
				Dalpha[5] = bv * (ezt * GND[2] + eyt * GND[3]);
				Dalpha[6] = bv * (ext * GND[3] + ezt * GND[1]);

				for (int i = 1; i < NDLOC; i++)
					for (int k = 1; k < NDS; k++)
						DISLOCATIONF[i] = DISLOCATIONF[i]
								+ (BBTRANE[i][k] * Dalpha[k] * Det);

			}

	AssemDislocationForce(IEL);

}

void DISLOCATIONELEMENTFORCE(double XLOC[NODEN], double YLOC[NODEN],
		double ZLOC[NODEN], int i) {

	if ((Xglob[LOTOGO[8 * i + (1 - 1)]-1] <= 0.0)
			&& (Xglob[LOTOGO[8 * i + (2 - 1)]-1] <= 0.0)
			&& (Xglob[LOTOGO[8 * i + (3 - 1)]-1] <= 0.0)
			&& (Xglob[LOTOGO[8 * i + (4 - 1)]-1] <= 0.0)
			&& (Xglob[LOTOGO[8 * i + (5 - 1)]-1] <= 0.0)
			&& (Xglob[LOTOGO[8 * i + (6 - 1)]-1] <= 0.0)
			&& (Xglob[LOTOGO[8 * i + (7 - 1)]-1] <= 0.0)
			&& (Xglob[LOTOGO[8 * i + (8 - 1)]-1] <= 0.0)) {
		bool condition1, condition2;

		condition1 = false;
		condition2 = false;

		if ((Yglob[LOTOGO[8 * i + (1 - 1)]-1] < 0.0)
				|| (Yglob[LOTOGO[8 * i + (2 - 1)]-1] < 0.0)
				|| (Yglob[LOTOGO[8 * i + (3 - 1)]-1] < 0.0)
				|| (Yglob[LOTOGO[8 * i + (4 - 1)]-1] < 0.0)
				|| (Yglob[LOTOGO[8 * i + (5 - 1)]-1] < 0.0)
				|| (Yglob[LOTOGO[8 * i + (6 - 1)]-1] < 0.0)
				|| (Yglob[LOTOGO[8 * i + (7 - 1)]-1] < 0.0)
				|| (Yglob[LOTOGO[8 * i + (8 - 1)]-1] < 0.0))
			condition1 = true;

		if ((Yglob[LOTOGO[8 * i + (1 - 1)]-1] > 0.0)
				|| (Yglob[LOTOGO[8 * i + (2 - 1)]-1] > 0.0)
				|| (Yglob[LOTOGO[8 * i + (3 - 1)]-1] > 0.0)
				|| (Yglob[LOTOGO[8 * i + (4 - 1)]-1] > 0.0)
				|| (Yglob[LOTOGO[8 * i + (5 - 1)]-1] > 0.0)
				|| (Yglob[LOTOGO[8 * i + (6 - 1)]-1] > 0.0)
				|| (Yglob[LOTOGO[8 * i + (7 - 1)]-1] > 0.0)
				|| (Yglob[LOTOGO[8 * i + (8 - 1)]-1] > 0.0))
			condition2 = true;

		if (condition1 && condition2) {

			double et[ND];

		
			et[0] = 0.0;
			et[1] = 0.0;
			et[2] = 0.0;
			et[3] = 1.0;

			for (int j = 1; j < NODEN; j++) {
				XLOC[j] = Xglob[LOTOGO[8 * i + (j - 1)]-1];
				YLOC[j] = Yglob[LOTOGO[8 * i + (j - 1)]-1];
				ZLOC[j] = Zglob[LOTOGO[8 * i + (j - 1)]-1];
			}

			DISLOCATIONFORCE(XLOC, YLOC, ZLOC, i, bv, et);
		}
	}
}

void STIFF(double XLOC[NODEN], double YLOC[NODEN], double ZLOC[NODEN]) {

	double A[NDS][10], B[10][10], C[10][NDLOC];
	double BB[NDS][NDLOC], BBTRAN[NDLOC][NDS];

	double SND[ND][NODEN];

	double BBTRANE[NDLOC][NDS], AB[NDS][10];

	double SHI, NEW, KSHI;
	double DIR1[NODEN], DIR2[NODEN], DIR3[NODEN], COR[NODEN][ND];

	DIR1[1] = DIR2[4] = DIR3[3] = -1.0;
	DIR1[2] = DIR2[8] = DIR3[4] = -1.0;
	DIR1[3] = DIR2[3] = DIR3[7] = 1.0;
	DIR1[4] = DIR2[7] = DIR3[8] = 1.0;
	DIR1[5] = DIR2[1] = DIR3[2] = -1.0;
	DIR1[6] = DIR2[5] = DIR3[1] = -1.0;
	DIR1[7] = DIR2[2] = DIR3[5] = 1.0;
	DIR1[8] = DIR2[6] = DIR3[6] = 1.0;

	for (int i = 0; i < NDS; i++)
		for (int j = 0; j < 10; j++)
			A[i][j] = 0.0;

	A[1][1] = A[2][5] = A[3][9] = 1.0;
	A[4][2] = A[4][4] = 1.0;
	A[5][6] = A[5][8] = 1.0;
	A[6][3] = A[6][7] = 1.0;

	for (int i = 0; i < NDLOC; i++)
		for (int j = 0; j < NDLOC; j++)
			AKLOC[i][j] = 0.0;

	for (int III = 1; III < 3; III++)        		// GAUSS INTEGRATION LOOP
		for (int JJJ = 1; JJJ < 3; JJJ++)
			for (int PPP = 1; PPP < 3; PPP++) {
				if (III == 1)
					SHI = -0.577350269189626;
				if (III == 2)
					SHI = 0.577350269189626;

				if (JJJ == 1)
					NEW = -0.577350269189626;
				if (JJJ == 2)
					NEW = 0.577350269189626;

				if (PPP == 1)
					KSHI = -0.577350269189626;
				if (PPP == 2)
					KSHI = 0.577350269189626;

				for (int i = 0; i < ND; i++)
					for (int j = 0; j < NODEN; j++)
						SND[i][j] = 0.0;

				for (int j = 1; j < NODEN; j++) {
					//Order changed to match Abaqus numbering
					SND[2][j] = (1.0 / 8.0) * DIR1[j] * (1.0 + NEW * DIR2[j])
							* (1.0 + KSHI * DIR3[j]);
					SND[1][j] = (1.0 / 8.0) * DIR2[j] * (1.0 + SHI * DIR1[j])
							* (1.0 + KSHI * DIR3[j]);
					SND[3][j] = (1.0 / 8.0) * DIR3[j] * (1.0 + NEW * DIR2[j])
							* (1.0 + SHI * DIR1[j]);
				}

				for (int i = 0; i < NODEN; i++)
					for (int j = 0; j < ND; j++)
						COR[i][j] = 0.0;

				for (int i = 1; i < NODEN; i++) {
					COR[i][1] = XLOC[i];
					COR[i][2] = YLOC[i];
					COR[i][3] = ZLOC[i];
				}

				for (int i = 0; i < ND; i++)
					for (int j = 0; j < ND; j++) {
						Ajac[i][j] = 0.0;
						AjacInv[i][j] = 0.0;
					}

				for (int i = 1; i < ND; i++)
					for (int j = 1; j < ND; j++)
						for (int k = 1; k < NODEN; k++)
							Ajac[i][j] = Ajac[i][j] + SND[i][k] * COR[k][j];

				JacInv();

				for (int i = 0; i < 10; i++)
					for (int j = 0; j < 10; j++)
						B[i][j] = 0.0;

				for (int j = 1; j < ND; j++)
					for (int l = 1; l < ND; l++) {
						B[j][l] = AjacInv[j][l];
						B[j + 3][l + 3] = AjacInv[j][l];
						B[j + 6][l + 6] = AjacInv[j][l];
					}

				for (int i = 0; i < 10; i++)
					for (int j = 0; j < NDLOC; j++)
						C[i][j] = 0.0;

				for (int i = 1; i < 4; i++)
					for (int j = 1; j < 24; j = j + 3)
						C[i][j] = SND[i][(j + 2) / 3];

				for (int i = 4; i < 7; i++)
					for (int j = 2; j < 24; j = j + 3)
						C[i][j] = SND[i - 3][(j + 1) / 3];

				for (int i = 7; i < 10; i++)
					for (int j = 3; j < 25; j = j + 3)
						C[i][j] = SND[i - 6][j / 3];

				for (int i = 0; i < NDS; i++)
					for (int j = 0; j < 10; j++)
						AB[i][j] = 0.0;

				for (int i = 1; i < NDS; i++)
					for (int j = 1; j < 10; j++)
						for (int k = 1; k < 10; k++)
							AB[i][j] = AB[i][j] + A[i][k] * B[k][j];

				for (int i = 0; i < NDS; i++)
					for (int j = 0; j < NDLOC; j++)
						BB[i][j] = 0.0;

				for (int i = 1; i < NDS; i++)
					for (int j = 1; j < NDLOC; j++)
						for (int k = 1; k < 10; k++)
							BB[i][j] = BB[i][j] + AB[i][k] * C[k][j];

				for (int i = 0; i < NDLOC; i++)
					for (int j = 0; j < NDS; j++) {
						BBTRAN[i][j] = 0.0;
						BBTRANE[i][j] = 0.0;
					}

				for (int i = 1; i < NDS; i++)
					for (int j = 1; j < NDLOC; j++)
						BBTRAN[j][i] = BB[i][j];

				for (int i = 1; i < NDLOC; i++)
					for (int j = 1; j < NDS; j++)
						for (int k = 1; k < NDS; k++)
							BBTRANE[i][j] = BBTRANE[i][j]
									+ BBTRAN[i][k] * EE3D[k][j];

				for (int i = 1; i < NDLOC; i++)
					for (int j = 1; j < NDLOC; j++)
						for (int k = 1; k < NDS; k++)
							AKLOC[i][j] = AKLOC[i][j]
									+ (BBTRANE[i][k] * BB[k][j] * Det);
			}
}


void STRESS(double XLOC[NODEN], double YLOC[NODEN], double ZLOC[NODEN],
		int IEL) {

	double A[NDS][10], B[10][10], C[10][NDLOC];
	double BB[NDS][NDLOC], BBTRAN[NDLOC][NDS];
	double SND[ND][NODEN];

	double BBTRANE[NDLOC][NDS], AB[NDS][10];

	double SHI, NEW, KSHI;
	double DIR1[NODEN], DIR2[NODEN], DIR3[NODEN], COR[NODEN][ND];
	double AvgStr1, AvgStr2, AvgStr3;
	
	bool condition1, condition2;

	DIR1[1] = DIR2[4] = DIR3[3] = -1.0;
	DIR1[2] = DIR2[8] = DIR3[4] = -1.0;
	DIR1[3] = DIR2[3] = DIR3[7] = 1.0;
	DIR1[4] = DIR2[7] = DIR3[8] = 1.0;
	DIR1[5] = DIR2[1] = DIR3[2] = -1.0;
	DIR1[6] = DIR2[5] = DIR3[1] = -1.0;
	DIR1[7] = DIR2[2] = DIR3[5] = 1.0;
	DIR1[8] = DIR2[6] = DIR3[6] = 1.0;

	for (int i = 0; i < NDS; i++)
		for (int j = 0; j < 10; j++)
			A[i][j] = 0.0;

	A[1][1] = A[2][5] = A[3][9] = 1.0;
	A[4][2] = A[4][4] = 1.0;
	A[5][6] = A[5][8] = 1.0;
	A[6][3] = A[6][7] = 1.0;

	/*for (int i = 0; i < NDLOC; i++)
		for (int j = 0; j < NDLOC; j++)
			AKLOC[i][j] = 0.0;*/

	AvgStr1 = 0.0;
	AvgStr2 = 0.0;
	AvgStr3 = 0.0;

	//cout << "RESULTS FOR ELEMENT#: " << IEL << endl;

	for (int III = 1; III < 3; III++)        		// GAUSS INTEGRATION LOOP
		for (int JJJ = 1; JJJ < 3; JJJ++)
			for (int PPP = 1; PPP < 3; PPP++) {
				if (III == 1)
					SHI = -0.577350269189626;
				if (III == 2)
					SHI = 0.577350269189626;

				if (JJJ == 1)
					NEW = -0.577350269189626;
				if (JJJ == 2)
					NEW = 0.577350269189626;

				if (PPP == 1)
					KSHI = -0.577350269189626;
				if (PPP == 2)
					KSHI = 0.577350269189626;

				for (int i = 0; i < ND; i++)
					for (int j = 0; j < NODEN; j++)
						SND[i][j] = 0.0;

				for (int j = 1; j < NODEN; j++) {
					//Order changed to match Abaqus numbering
					SND[2][j] = (1.0 / 8.0) * DIR1[j] * (1.0 + NEW * DIR2[j])
							* (1.0 + KSHI * DIR3[j]);
					SND[1][j] = (1.0 / 8.0) * DIR2[j] * (1.0 + SHI * DIR1[j])
							* (1.0 + KSHI * DIR3[j]);
					SND[3][j] = (1.0 / 8.0) * DIR3[j] * (1.0 + NEW * DIR2[j])
							* (1.0 + SHI * DIR1[j]);
				}

				for (int i = 0; i < NODEN; i++)
					for (int j = 0; j < ND; j++)
						COR[i][j] = 0.0;

				for (int i = 1; i < NODEN; i++) {
					COR[i][1] = XLOC[i];
					COR[i][2] = YLOC[i];
					COR[i][3] = ZLOC[i];
				}

				for (int i = 0; i < ND; i++)
					for (int j = 0; j < ND; j++) {
						Ajac[i][j] = 0.0;
						AjacInv[i][j] = 0.0;
					}

				for (int i = 1; i < ND; i++)
					for (int j = 1; j < ND; j++)
						for (int k = 1; k < NODEN; k++)
							Ajac[i][j] = Ajac[i][j] + SND[i][k] * COR[k][j];

				JacInv();

				for (int i = 0; i < 10; i++)
					for (int j = 0; j < 10; j++)
						B[i][j] = 0.0;

				for (int j = 1; j < ND; j++)
					for (int l = 1; l < ND; l++) {
						B[j][l] = AjacInv[j][l];
						B[j + 3][l + 3] = AjacInv[j][l];
						B[j + 6][l + 6] = AjacInv[j][l];
					}

				for (int i = 0; i < 10; i++)
					for (int j = 0; j < NDLOC; j++)
						C[i][j] = 0.0;

				for (int i = 1; i < 4; i++)
					for (int j = 1; j < 24; j = j + 3)
						C[i][j] = SND[i][(j + 2) / 3];

				for (int i = 4; i < 7; i++)
					for (int j = 2; j < 24; j = j + 3)
						C[i][j] = SND[i - 3][(j + 1) / 3];

				for (int i = 7; i < 10; i++)
					for (int j = 3; j < 25; j = j + 3)
						C[i][j] = SND[i - 6][j / 3];

				for (int i = 0; i < NDS; i++)
					for (int j = 0; j < 10; j++)
						AB[i][j] = 0.0;

				for (int i = 1; i < NDS; i++)
					for (int j = 1; j < 10; j++)
						for (int k = 1; k < 10; k++)
							AB[i][j] = AB[i][j] + A[i][k] * B[k][j];

				for (int i = 0; i < NDS; i++)
					for (int j = 0; j < NDLOC; j++)
						BB[i][j] = 0.0;

				for (int i = 1; i < NDS; i++)
					for (int j = 1; j < NDLOC; j++)
						for (int k = 1; k < 10; k++)
							BB[i][j] = BB[i][j] + AB[i][k] * C[k][j];

				double DIS[NODEN][ND], DP[NDLOC], BBDP[NDS], ENDP[NDLOC];
				int iicon;

				for (int j = 1; j < NODEN; j++)
					for (int jjk = 1; jjk < ND; jjk++) {
					    iicon = JJglob[3 * LOTOGO[8 * IEL + (j - 1)] + (jjk - 1)]; 
						if (iicon < 0)
							DIS[j][jjk] = 0.0;  // ??? constraint value?
						else
							DIS[j][jjk] = Dglob[iicon];
					}

				for (int i = 1; i < NODEN; i++) {
					DP[(3 * i - 2)] = DIS[i][1];
					DP[(3 * i - 1)] = DIS[i][2];
					DP[3 * i] = DIS[i][3];
				}

				for (int i = 1; i < NODEN; i++) {
					ENDP[(3 * i - 2)] = ENFRDISPglob[3
							* LOTOGO[8 * IEL + (i - 1)] + (1 - 1)];
					ENDP[(3 * i - 1)] = ENFRDISPglob[3
							* LOTOGO[8 * IEL + (i - 1)] + (2 - 1)];
					ENDP[3 * i] = ENFRDISPglob[3
					        * LOTOGO[8 * IEL + (i - 1)] + (3 - 1)];
				}

				for (int i = 1; i < NDS; i++) {
					STR[i] = 0.0;
					BBDP[i] = 0.0;
				}

                 int i = IEL;
	        

        if ((Xglob[LOTOGO[8 * i + (1 - 1)]-1] <= 0.0)
			&& (Xglob[LOTOGO[8 * i + (2 - 1)]-1] <= 0.0)
			&& (Xglob[LOTOGO[8 * i + (3 - 1)]-1] <= 0.0)
			&& (Xglob[LOTOGO[8 * i + (4 - 1)]-1] <= 0.0)
			&& (Xglob[LOTOGO[8 * i + (5 - 1)]-1] <= 0.0)
			&& (Xglob[LOTOGO[8 * i + (6 - 1)]-1] <= 0.0)
			&& (Xglob[LOTOGO[8 * i + (7 - 1)]-1] <= 0.0)
			&& (Xglob[LOTOGO[8 * i + (8 - 1)]-1] <= 0.0)) {

		condition1 = false;
		condition2 = false;

		if ((Yglob[LOTOGO[8 * i + (1 - 1)]-1] < 0.0)
				|| (Yglob[LOTOGO[8 * i + (2 - 1)]-1] < 0.0)
				|| (Yglob[LOTOGO[8 * i + (3 - 1)]-1] < 0.0)
				|| (Yglob[LOTOGO[8 * i + (4 - 1)]-1] < 0.0)
				|| (Yglob[LOTOGO[8 * i + (5 - 1)]-1] < 0.0)
				|| (Yglob[LOTOGO[8 * i + (6 - 1)]-1] < 0.0)
				|| (Yglob[LOTOGO[8 * i + (7 - 1)]-1] < 0.0)
				|| (Yglob[LOTOGO[8 * i + (8 - 1)]-1] < 0.0))
			condition1 = true;

		if ((Yglob[LOTOGO[8 * i + (1 - 1)]-1] > 0.0)
				|| (Yglob[LOTOGO[8 * i + (2 - 1)]-1] > 0.0)
				|| (Yglob[LOTOGO[8 * i + (3 - 1)]-1] > 0.0)
				|| (Yglob[LOTOGO[8 * i + (4 - 1)]-1] > 0.0)
				|| (Yglob[LOTOGO[8 * i + (5 - 1)]-1] > 0.0)
				|| (Yglob[LOTOGO[8 * i + (6 - 1)]-1] > 0.0)
				|| (Yglob[LOTOGO[8 * i + (7 - 1)]-1] > 0.0)
				|| (Yglob[LOTOGO[8 * i + (8 - 1)]-1] > 0.0))
			condition2 = true;

		if (condition1 && condition2) {
                        
			for (int j = 1; j < NODEN; j++){
                if (Yglob[LOTOGO[8 * i + (j - 1)]-1] > 0)
					DP[3 * j] = DP[3 * j] - 0.5 * bv;
                else
                    DP[3 * j] = DP[3 * j] + 0.5 * bv;
				}
                        
        }
	    }


				for (int i = 1; i < NDS; i++)
					for (int k = 1; k < NDLOC; k++)
						BBDP[i] = BBDP[i] + BB[i][k] * (DP[k] + ENDP[k]);

				BBDP[4] = 2.0 * BBDP[4];
                                BBDP[5] = 2.0 * BBDP[5];
                                BBDP[6] = 2.0 * BBDP[6];

				for (int i = 1; i < NDS; i++)
					for (int k = 1; k < NDS; k++)
						STR[i] = STR[i] + EE3D[i][k] * BBDP[k];

				AvgStr1 = AvgStr1 + 0.125 * STR[1];
				AvgStr2 = AvgStr2 + 0.125 * STR[5];
				AvgStr3 = AvgStr3 + 0.125 * STR[6];
			}

    sigma[0] = AvgStr1;
    sigma[1] = AvgStr2;
    sigma[2] = AvgStr3;

}


void create_mesh() {
	/* ******************** BEGIN FEM ELEMENT MESH CREATION ******************** */

	int indx, indy, indz, ii;
	double valx, valy, valz, w;
	double distx, disty, distz;
	double displacementx, displacementy, displacementz;
	double newdist_fem[3];
	double newdist_atom[3];
	int count_elems_fem_x, count_elems_atom_x;
	int count_elems_fem_y, count_elems_atom_y;
	nelem_full_big_box_x = (fem_elems_round[0]/2)*2 + num_space_x;
	nelem_full_big_box_y = (fem_elems_round[1]/2)*2 + num_space_y;
	nelem_full_big_box_z = (Lz+0.02)/(2*dist[2]);
	bool skip_elem_x = 0;
	bool skip_elem_y = 0;
	newdist_fem[0] = (fem_big_box[0]-(0.5*sys_width[0]-0.05))/(fem_elems_round[0]/2);
	newdist_fem[1] = (fem_big_box[1]-(0.5*sys_width[1]-0.05))/(fem_elems_round[1]/2);
	newdist_fem[2] = dist[2];
	newdist_atom[0] = sys_reduced[0]/num_space_x;
	newdist_atom[1] = sys_reduced[1]/num_space_y;
	newdist_atom[2] = dist[2];
	int nodecount = 0;
	int elemcount = 0;
	int innersurfacecount = 0;
	int outersurfacecount = 0;
	int backsurfacecount = 0;
	cout << "newdist_fem is (" << newdist_fem[0] << "," << newdist_fem[1] << "," << newdist_fem[2] << ") ; newdist_atom is (" << newdist_atom[0] << "," << newdist_atom[1] << "," << newdist_atom[2] << ")" << endl;
	distx = newdist_fem[0];
	distz = newdist_fem[2];
	displacementx = 0.0;
	displacementy = 0.0;
	displacementz = 0.0;
	count_elems_fem_x = 0;
	count_elems_atom_x = 0;

	// correction variable, which will be modified to parse nodes with an "x" that falls in the central hole
	double shift_y = nelem_full_big_box_y;
	for (indx = 0; indx <= nelem_full_big_box_x; indx++) {
		skip_elem_x = 0;
		valx = -fem_big_box[0] + (indx-count_elems_fem_x-count_elems_atom_x)*distx + count_elems_fem_x*newdist_fem[0] + count_elems_atom_x*newdist_atom[0];
		if (fabs(fabs(valx)-0.5*sys_reduced[0]) < 0.02) {
			if (valx < 0) {
				distx = newdist_atom[0];
				count_elems_fem_x = indx;
				skip_elem_x = 1;
			}
			else {
				distx = newdist_fem[0];
				count_elems_atom_x = indx-count_elems_fem_x;
			}
		}
		if (fabs(valx) < (0.5*sys_reduced[0]-0.02)) {
			skip_elem_x = 1;
			// initially correct value of "shift_y" here
			shift_y = nelem_full_big_box_y - num_space_y + 1;
		}
		// indy less or equal in order to reach valy = -fem_big_box[1]
		disty = newdist_fem[1];
		count_elems_fem_y = 0;
		count_elems_atom_y = 0;
		for (indy = 0; indy <= nelem_full_big_box_y; indy++) {
			skip_elem_y = 0;
			valy = -fem_big_box[1] + (indy-count_elems_fem_y-count_elems_atom_y)*disty + count_elems_fem_y*newdist_fem[1] + count_elems_atom_y*newdist_atom[1];
			if (fabs(fabs(valy)-0.5*sys_reduced[1]) < 0.02) {
				if (valy < 0) {
					disty = newdist_atom[1];
					count_elems_fem_y = indy;
					skip_elem_y = 1;
				}
				else {
					disty = newdist_fem[1];
					count_elems_atom_y = indy-count_elems_fem_y;
				}
			}
			if (fabs(valy) < (0.5*sys_reduced[1]-0.02)) {
				skip_elem_y = 1;
			}
			// correct "shift_y" for low "x" values
			if (fabs(valx+0.5*sys_reduced[0]) < 0.02 && (valy-0.5*sys_reduced[1]+0.02) > 0) {
				shift_y = nelem_full_big_box_y - num_space_y + 1;
			}
			// un-correct "shift_y" for high "x" values
			if (fabs(valx-(0.5*sys_reduced[0]-newdist_atom[0])) < 0.02 && (valy-0.5*sys_reduced[1]+0.02) > 0) {
				shift_y = nelem_full_big_box_y;
			}
			if (fabs(valx) < (0.5*sys_reduced[0]-0.02) && fabs(valy) < (0.5*sys_reduced[1]-0.02)) continue;
			for (indz = 0; indz <= nelem_full_big_box_z; indz++) {
				valz = indz*2*dist[2];
				nodecount++;
				Xglob.push_back(valx);
				Yglob.push_back(valy);
				Zglob.push_back(valz);
				if (indz == 0) {
					backsurfacecount++;
				}
				if (indx < nelem_full_big_box_x && indy < nelem_full_big_box_y && indz < nelem_full_big_box_z && !(skip_elem_x && skip_elem_y)) {
					elemcount++;
					LOTOGO.push_back(nodecount);
					LOTOGO.push_back(nodecount + nelem_full_big_box_z + 1);
					LOTOGO.push_back(nodecount + nelem_full_big_box_z + 2);
					LOTOGO.push_back(nodecount + 1);
					LOTOGO.push_back(nodecount + (shift_y+1)*(nelem_full_big_box_z+1));
					LOTOGO.push_back(nodecount + (shift_y+1)*(nelem_full_big_box_z+1) 
                            					   + nelem_full_big_box_z + 1);
					LOTOGO.push_back(nodecount + (shift_y+1)*(nelem_full_big_box_z+1)
                            					   + nelem_full_big_box_z + 2);
					LOTOGO.push_back(nodecount + (shift_y+1)*(nelem_full_big_box_z+1) + 1);
				}
				if ((fabs(fabs(valx)-0.5*sys_reduced[0]) < 0.02 && fabs(valy) < 0.5*sys_reduced[1]+0.02) || (fabs(valx) < 0.5*sys_reduced[0]+0.02 && fabs(fabs(valy)-0.5*sys_reduced[1]) < 0.02)) {
                    nID[innersurfacecount][0] = nodecount;
                    nID[innersurfacecount][1] = valx;
                    nID[innersurfacecount][2] = valy;
                    nID[innersurfacecount][3] = valz;
					innersurfacecount++;
        			ncstr_i.push_back(nodecount);
        				}
				if ((fabs(fabs(valx)-fem_big_box[0]) < 0.02) || (fabs(fabs(valy)-fem_big_box[1]) < 0.02)) {
					outersurfacecount++;
        			ncstr_o.push_back(nodecount);
        				}
			   }
		}
	}

	NLAG = backsurfacecount-(innersurfacecount+outersurfacecount)/(nelem_full_big_box_z+1);
	NNODE = nodecount;
	NEL = elemcount;
	NCF = 0;
	NCONT = innersurfacecount + outersurfacecount;
	NCONSNODE = innersurfacecount;
	NENFD = outersurfacecount;
	node_dis.resize(3*NNODE, 0.0);
		
	JJglob.resize((3 * (NNODE+1)), 1);
	JJglob[0] = -1;
    JJglob[1] = -1;
    JJglob[2] = -1;
	ENFRDISPglob.resize((3 * (NNODE+1)), 0.0);
	int inode;
	for (ii = 0; ii < NCONSNODE; ++ii) {
        // set conversion matrix between nodal dof and solution vector
        inode = ncstr_i[ii];
        JJglob[3*inode  ] = -1;
        JJglob[3*inode+1] = -1;
        JJglob[3*inode+2] = -1;
        }
	for (ii = 0; ii < NENFD; ++ii) {
        // set conversion matrix between nodal dof and solution vector
        inode = ncstr_o[ii];
        JJglob[3*inode  ] = -1;
        JJglob[3*inode+1] = -1;
        JJglob[3*inode+2] = -1;
        }
	NDF = EqNum(&JJglob[0], NNODE);
	Fglob.resize((NDF), 0.0);
	xdof = NDF + 3 * NLAG;
	
	// select nodal_pairs_along_dislocation_line
    int counter = 0;
	for (ii = 0; ii < NNODE; ++ii){
		if (ii % (nelem_full_big_box_z+1) != 0) continue;
		if (!((fabs(fabs(Xglob[ii])-0.5*sys_reduced[0]) < 0.02 && 
    		       fabs(Yglob[ii]) < 0.5*sys_reduced[1]+0.02) || 
    		      (fabs(Xglob[ii]) < 0.5*sys_reduced[0]+0.02 &&
    		       fabs(fabs(Yglob[ii])-0.5*sys_reduced[1]) < 0.02)) &&
		    !((fabs(fabs(Xglob[ii])-fem_big_box[0]) < 0.02) ||
		      (fabs(fabs(Yglob[ii])-fem_big_box[1]) < 0.02))) {
			ContraintNodes.push_back(ii + 1);
			ContraintNodes.push_back(ii + nelem_full_big_box_z + 1);
			++counter;
        }
    }
    if (NLAG != counter) {
        cout << "WARNING!!!!!!!" << endl;
        cout << "NLAG: " << NLAG << " Nodal pairs: " << counter << endl;
    }
	
}

void lagrange_constraint() {

	int icon, jcon, iicon, jjcon, ind;
	
	for (int i = 0; i < NLAG; i++){
        icon = ContraintNodes[2*i    ];
        jcon = ContraintNodes[2*i + 1];
        for (int j = 0; j < 3; j++) {
        
            iicon = JJglob[3*icon + j];
            jjcon = JJglob[3*jcon + j];
    
            if ((iicon < 0) || (jjcon < 0)) {
                cout << "Error in Lagrange constraint !" << endl;
                cout << "LOOP: " << i << " in: " << j << " NLAG: " << NLAG << endl;
                cout << "Node pair: " << icon << ", " << jcon << endl;
                cout << "EQ nrs: " << iicon << ", " << jjcon << endl;
                exit(EXIT_FAILURE);
            }
    
            MAT[iicon][NDF + 3*i + j] = 1;
            MAT[jjcon][NDF + 3*i + j] = -1;
        }
	}

	for (ind = NDF; ind < xdof; ind++) {
        auto it = MAT[ind].find(ind);

        if (it == MAT[ind].end())
            MAT[ind][ind] = 0;
        else {
            cout << "Error in Lagrange matrix addition!" << endl;
            exit(EXIT_FAILURE);
      }
	}
}

void Assem(int Iel) {
		/* Insert element stiffness matrix into system stiffness
		matrix MAT */
	int i, k, icon, kcon, iicon, kkcon, ii, kk, iii, kkk, ilp;

	MAT.resize(NDF + 3 * NLAG);

	for (i = 1; i < NODEN; i++)
		for (k = 1; k < NODEN; k++) {
			icon = LOTOGO[8 * Iel + (i - 1)];
			kcon = LOTOGO[8 * Iel + (k - 1)];
			for (ii = 1; ii < ND; ii++)
				for (kk = 1; kk < ND; kk++) {
					iicon = JJglob[(3 * icon) + (ii - 1)];
					kkcon = JJglob[(3 * kcon) + (kk - 1)];
					iii = 3 * (i - 1) + ii;
					kkk = 3 * (k - 1) + kk;
					if ((iicon >= 0) && (kkcon >= 0)) {
						if (iicon <= kkcon) {
                            // only store upper triangle of symmetric matrix
							auto it = MAT[iicon].find(kkcon);

							if (it != MAT[iicon].end()) {
								MAT[iicon][kkcon] += AKLOC[iii][kkk];
							} else {
								MAT[iicon][kkcon] = AKLOC[iii][kkk];
							}

						}
					}
				}
		}
}

void skylineTOcompressedarray() {
    /* Transform system stiffness matrix MAT from skyline format into
    Compressed Sparse Row (CSR) format with data in NEWAglob, columns in Iglob and 
    rows in Jglob*/
	int colht, count, Icount, IGcount;
	double val, Ival, Jval;

	IGcount = count = 0;

	int j = 0;
	double value;
	Iglob.push_back(0);
	for (int i = 0; i < xdof; ++i) {
		for (auto it = MAT[i].begin(); it != MAT[i].end(); ++it) {
			val = it->second;
			if ((val==0.0) || (fabs(val) > 1.e-9)) {
			    NEWAglob.push_back(val);
			    Jglob.push_back(it->first);
			    ++j;
			}
		}
		Iglob.push_back(j);
	}
    int size_Aglob = NEWAglob.size(); 
	cout << size_Aglob << "  " << Iglob.size() << "   " << Jglob.size() << endl;
    cout << " n = " << NDF + 3 * NLAG << " nmax = " << size_Aglob << endl;
}

void init_matrix() {
    /* Assemble stiffness matrix
    */
    double XLOC[NODEN], YLOC[NODEN], ZLOC[NODEN];
	cout << "LOOPING OVER ELEMENTS " << endl;

	for (int IEL = 0; IEL < NEL; IEL++) {
		for (int j = 0; j < 8; j++) {
			XLOC[j+1] = Xglob[LOTOGO[8 * IEL + j] - 1];
			YLOC[j+1] = Yglob[LOTOGO[8 * IEL + j] - 1];
			ZLOC[j+1] = Zglob[LOTOGO[8 * IEL + j] - 1];
		}

		STIFF(XLOC, YLOC, ZLOC);
		Assem(IEL);
	}

	cout << "COMPUTING COMPRESSED ARRAY " << endl;
	lagrange_constraint();
	skylineTOcompressedarray();
}

void update_fglob() {
    /* Update global force vector
    */
    double XLOC[NODEN], YLOC[NODEN], ZLOC[NODEN];

	for (int IEL = 0; IEL < NEL; IEL++) {
		for (int j = 0; j < 8; j++) {
			XLOC[j+1] = Xglob[LOTOGO[8 * IEL + j] - 1];
			YLOC[j+1] = Yglob[LOTOGO[8 * IEL + j] - 1];
			ZLOC[j+1] = Zglob[LOTOGO[8 * IEL + j] - 1];
		}

		ENFORCEDFORCE(XLOC, YLOC, ZLOC, IEL);
		Assem_force(IEL);
		DISLOCATIONELEMENTFORCE(XLOC, YLOC, ZLOC, IEL);
	}
}

void create_volterra_dis() {
    /* create BC for Volterra dislocation on constraint nodes
    
    Parameters:

    */
    int i, inode;
    double yc;
    
    // outer boundary
    for (i = 0; i < NENFD; i++) {
		inode = ncstr_o[i];
		ENFRDISPglob[3*inode    ] = 0.0;
		ENFRDISPglob[3*inode + 1] = 0.0;
		ENFRDISPglob[3*inode + 2] = bv * atan2(Yglob[inode-1], Xglob[inode-1])/(2.0*PI);
	}	
	//inner boundary
	for (i = 0; i < NCONSNODE; i++) {
	    inode = ncstr_i[i];
		ENFRDISPglob[3*inode    ] = 0.0;
		ENFRDISPglob[3*inode + 1] = 0.0;
	    ENFRDISPglob[3*inode + 2] = bv * atan2(Yglob[inode-1], Xglob[inode-1])/(2.0*PI);
	}
}

void apply_e23_outer(double e23) {
    /* apply given yz-strain component as displacement bc to outer nodes
    
    Parameters:
    e23 : double
        yz-shear strain
    */
    int i, inode;
    double yc;
    
    for (i = 0; i < NENFD; i++) {
		inode = ncstr_o[i];
		yc = Yglob[inode-1];
		//ENFRDISPglob[3*inode    ] += 0.0;
		//ENFRDISPglob[3*inode + 1] += 0.0;
		ENFRDISPglob[3*inode + 2] += yc * e23;
	}
}
