






// void InsertOneDimMaterial(TPZGeoMesh *gmesh)
// {
// 	// Inserir elmentos 1D fmatLambda and fmatLambdaBCs
// 	int64_t nel = gmesh->NElements();

// 	for (int64_t el = 0; el < nel; el++)
// 	{
// 		TPZGeoEl *gel = gmesh->Element(el);
// 		if (gel->HasSubElement() && f_allrefine){continue;}
// 		if (gel->Dimension() != gmesh->Dimension()){continue;}
// 		int nsides = gel->NSides();
// 		for (int is = 0; is < nsides; is++)
// 		{
// 			if (gel->SideDimension(is) != gmesh->Dimension() - 1){continue;}
// 			TPZGeoElSide gelside(gel, is);
// 			TPZGeoElSide neighbour = gelside.Neighbour();
// 			if (neighbour == gelside && f_allrefine == false){continue;}

// 			while (neighbour != gelside)
// 			{
// 				if (neighbour.Element()->Dimension() == gmesh->Dimension() - 1)
// 				{
// 					int neigh_matID = neighbour.Element()->MaterialId();
// 					if (neigh_matID == fmatBCbott){
// 						TPZGeoElBC(gelside, fmatLambdaBC_bott);
// 					}
// 					else if (neigh_matID == fmatBCtop){
// 						TPZGeoElBC(gelside, fmatLambdaBC_top);
// 					}
// 					else if (neigh_matID == fmatBCleft){
// 						TPZGeoElBC(gelside, fmatLambdaBC_left);
// 					}
// 					else if (neigh_matID == fmatBCright){
// 						TPZGeoElBC(gelside, fmatLambdaBC_right);
// 					}
// 					break;
// 				}
// 				if (neighbour.Element()->HasSubElement()){
// 					break;
// 				}
// 				neighbour = neighbour.Neighbour();
// 			}

// 			if (neighbour == gelside)
// 			{
// 				TPZGeoElBC(gelside, fmatLambda);
// 			}
// 		}
// 	}
// }