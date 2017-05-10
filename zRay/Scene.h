/*
 * zRay ~ A Monte Carlo Path Tracer
 * Zephyr Mann ~ 2013
 * SCAD ~ Studio 2
 * zephmann@gmail.com
 *
 * Scene.h
 */

#ifndef SCENE_H
#define SCENE_H

class Scene{
public:
	Scene();

	bool loadScene(string filename, string &outfile, Camera &c, int &w, int &h, bool &fast, int &samps, int &aadivs, double &noiseThreshold);
	Pixel Render(Ray r, int E=1) const;

	void Test() {
		//cout << numShapes << endl;
		//shapes[0]->Test();
		cout << mats[1].diffTex << endl;
	}

private:
	Shape **shapes;
	int numShapes;
	Light *lights;
	int numLights;
	Material *mats;
	int numMats;
	int *emitters;
	int numEmit;
	Texture *tex;
	int numTex;

	int maxdepth;
	int mindepth;
	int maxdiff;
	Pixel bgcolor;
};

Scene::Scene() {
}

Pixel Scene::Render(const Ray r, int E) const {
	Pixel p;

	double t = INFINITY;
	int i = -1;
	double tu,tv,u,v;
	int ttri,tri,bboxPlane;
	double tmin, tmax; //bounding box intersection distances

	// Find Intersection
	for (int ti = 0; ti < numShapes; ti++)
	{
		if( (r.depth > 0 || shapes[ti]->primary) &&
			(!r.direct || shapes[ti]->shadows) )
		{

			double tt,ttmin,ttmax;

			if(shapes[ti]->box.Intersect(r,ttmin,ttmax,bboxPlane) && shapes[ti]->Intersect(r,tt,tu,tv,ttri,ti==r.index,r.tri))
			{
				if(tt < t)
				{
					t = tt;
					i = ti;
					u = tu;
					v = tv;
					tri = ttri;
				}

				tmin = ttmin;
				tmax = ttmax;
			}
		}
	}

	// Shading routine
	if (i != -1)
	{
		bool isLight = false;
		for(int lig = 0; lig < numEmit; lig++)
		{
			if(i == emitters[lig])
				isLight = true;
		}

		Point interP = r(t);
		Normal n = shapes[i]->GetNormal(interP,tri,u,v);

		Material m = mats[shapes[i]->matId];

		double texU, texV;
		shapes[i]->GetUV(tri,u,v,texU,texV);

		Pixel ambi(m.amr, m.amg, m.amb);
		if(m.ambiTex != -1)
			ambi = tex[m.ambiTex].getPixel(texU, texV);

		if(isLight)
		{
			Pixel emit(m.er, m.eg, m.eb);
			if(m.emitTex != -1)
				emit = tex[m.emitTex].getPixel(texU, texV);

			p = ambi * emit;
		}
		else
			p = ambi;


		Pixel color(m.r, m.g, m.b);
		if(m.diffTex != -1)
			color = tex[m.diffTex].getPixel(texU, texV);

		//double falloff = 1. / pow(t,2);
		//p = p * falloff;

		double russR = color.r * m.diff;
		double russG = color.g * m.diff;
		double russB = color.b * m.diff;


		double rr = (russR > russG && russR > russB) ? russR : (russG > russB) ? russG : russB;

		if(r.depth < mindepth || (rr && r.depth < maxdepth && rand01() < rr))
		{

			bool inside = false;

			double opacity = m.opr;
			if(m.opacTex != -1)
				opacity = tex[m.opacTex].getGray(texU, texV);

			if(i == -5)
			{
				double ttu = fmod(texU*50,1.);
				double ttv = fmod(texV*50,1.);

				if (!((ttu > 0.5 && ttv > 0.5) || (ttu < 0.5 && ttv < 0.5)) ){
					color.r = 1. - color.r;
					color.g = 1. - color.g;
					color.b = 1. - color.b;
				}
			}

			if(Dot(n, r.d) > 0 )
			{
				n = n * -1;
				inside = true;
			}


			//offset normal
			//if(m.bumpTex != -1 && m.bumpDepth > 0.)
			if(m.bumpDepth > 0.)
			{
				Vector bump(0., 0., 0.);
				if(m.bumpTex != -1)
				{
					bump = Vector(tex[m.bumpTex].getPixel(texU, texV));
					bump = (bump * 2) - Vector(1,1,1);
				}
				else
				{
					double tNoise = pnoise(m.freqX * interP.x + m.offX,
						m.freqY * interP.y + m.offY,
						m.freqZ * interP.z + m.offZ);
					bump = Vector(tNoise, tNoise, tNoise);
				}


				bump = bump * PI_2 * m.bumpDepth;

				Vector w = Vector(n);
				Vector u = fabs(w.x) > 0.1 ? Vector(0,1,0) : Vector(1,0,0);
				u = Normalize(Cross(u,w));
				Vector v = Cross(w,u);

				Vector tn = w*cos(bump.x) + u*sin(bump.x);
				n = Normalize(Normal(tn*cos(bump.y) + v*sin(bump.y)));
			}

			// Diffuse
			if(m.diff > 0 && opacity > 0)
			{
				// DIRECT ILLUMINATION
				if(r.diffDepth < maxdiff)
				{
					//for(int j = 0; j < numEmit; j++)
					{
						int j = rand() % numEmit;

						Point emitSamp = shapes[emitters[j]]->GetRandomPoint();
						Vector emitDir = Normalize(emitSamp - interP);

						double emDot = Dot(emitDir, n);
						if(emDot > 0)
						{
							Pixel lp = Render(Ray(interP, emitDir, 0, INFINITY, r.depth+1, r.diffDepth+1, r.opacDepth, i, tri, true));
							//Pixel lp(mats[shapes[emitters[j]]->matId].er, mats[shapes[emitters[j]]->matId].eg, mats[shapes[emitters[j]]->matId].eb);

							lp = lp * emDot;
							lp = lp * color;
							lp = lp * m.diff;
							lp = lp * opacity;
							p = p + lp;
						}
					}
				}

				// INDIRECT ILLUMINATION
				double r1 = 2 * PI * rand01();
				double r2 = rand01();
				double r2s = sqrt(r2);

				Vector w = Vector(n);
				Vector u = fabs(w.x) > 0.1 ? Vector(0,1,0) : Vector(1,0,0);
				u = Normalize(Cross(u,w));
				Vector v = Cross(w,u);

				Vector psi = Normalize(u*cos(r1)*r2s + v*sin(r1)*r2s + w*sqrt(1-r2));

				Pixel rp = Render(Ray(interP, psi, 0, INFINITY, r.depth+1, r.diffDepth+1, r.opacDepth, i, tri), 1);

				rp = rp * color;
				rp = rp * m.diff;
				rp = rp * opacity;
				p = p + rp;
			}
			//End Diffuse

			//Specular
			double nOpac = r.opacDepth - (r.opacDepth * opacity);
			if(nOpac < EPSILON)
				nOpac = 0.;

			Vector reflV;
			double trefl = 1.;
			if(m.refl > 0 || m.refr)
			{
				trefl = Dot(r.d, -n);
				reflV = Vector(2. * trefl * n);
				reflV = Normalize(r.d + reflV);
			}

			Vector refrDir = r.d;
			double refrIndex = 1.;
			if(m.refr)
				refrIndex = inside ? 1. / m.refIndex : m.refIndex;

			double trefr = ( 1 - (refrIndex*refrIndex) * (1 - (trefl*trefl)) );

			refrDir = Normalize(Vector(r.d * refrIndex) + ((refrIndex * trefl - sqrtf(trefr)) * Vector(n)));

			double Re = 1.;
			double Tr = 1.;

			if(m.refl > 0 && nOpac > 0)
			{
				double a = refrIndex - 1;
				double b = refrIndex + 1;
				double R0 = (a*a) / (b*b);
				double c = 1 - (inside ? Dot(refrDir, -n) : trefl);

				//Re = R0 + (1-R0)*c*c*c*c*c;
				Re = R0 + (1-R0)*c*c;
				Tr = 1-Re;
			}

			// reflections
			Pixel rp;
			if(m.phong && m.refl > 0)
			{
				if(m.spec > 0)
				{
					double r1 = 2 * PI * rand01();
					double r2 = m.spec*rand01();
					double r2s = sqrt(r2);

					Vector w = reflV;
					Vector u = fabs(w.x) > 0.1 ? Vector(0,1,0) : Vector(1,0,0);
					u = Cross(u,w);
					Vector v = Cross(w,u);

					reflV = Normalize(u*cos(r1)*r2s + v*sin(r1)*r2s + w*sqrt(1-r2));
				}


				if(Dot(reflV, n) > 0)
				{
					Ray reflR(interP, reflV, 0, INFINITY, r.depth+1, r.diffDepth, r.opacDepth, i, tri);

					Pixel spec(m.sr, m.sg, m.sb);
					if(m.specTex != -1)
						spec = tex[m.specTex].getPixel(texU, texV);

					//Pixel rp = Render(reflR);
					rp = Render(reflR);
					rp = rp * spec;
					rp = rp * m.refl;

					if(m.refr && trefr >= 0)
						rp = rp * Re;

					p = p + rp;
				}
			}

			// transparency
			if(nOpac > 0 && (!m.refr || trefr >= 0))
			{
				Ray opacR(interP, refrDir, 0, INFINITY, r.depth+1, r.diffDepth, nOpac, i, tri);

				Pixel op = Render(opacR);

				op = op * r.opacDepth;
				op = op * Tr;

				p = p + op;

				if(0 && !m.refr)
				{
					p.r = max(op.r, rp.r);
					p.g = max(op.g, rp.g);
					p.b = max(op.b, rp.b);
				}
			}
			// End Specular

			// Color by Normals
			if(0) {
				p.r = n.x * 0.5 + 0.5;
				p.g = n.y * 0.5 + 0.5;
				p.b = n.z * 0.5 + 0.5;
			}

		} // if less than max depth

	} //if (index != 1)

	else
		p = bgcolor;

	return p;
}
// END RENDER //


bool Scene::loadScene(string filename, string &outfile, Camera &c, int &w, int &h, bool &fast, int &samps, int &aadivs, double &noiseThreshold) {

	//string filename = "C:\\images\\test.txt";
	ifstream sceneFile;
	sceneFile.open(filename.c_str(), ios::in);

	string word1,word2,word3;

	if(sceneFile.is_open())
	{
		double invRad = 1 / 180.; //used for conversion from degrees to radians

		sceneFile >> word1;  //scene name
		sceneFile >> word1;  //output file
		outfile = word1;

		// read width and height
		sceneFile >> word1 >> word2;
		w = atoi(word1.c_str());
		h = atoi(word2.c_str());

		// ray depth
		sceneFile >> word1 >> word2 >> word3;
		mindepth = atoi(word1.c_str());
		maxdepth = atoi(word2.c_str());
		maxdiff = atoi(word3.c_str());

		//
		sceneFile >> word1 >> word2 >> word3;
		fast = atoi(word1.c_str());
		aadivs = atoi(word2.c_str());
		samps = atoi(word3.c_str());

		sceneFile >> word1;
		noiseThreshold = atof(word1.c_str());

		// read in camera information
		sceneFile >> word1 >> word1;  //"Camera: Position"
		double x,y,z,aov;  // doubles to store position, rotation, and angle of view
		sceneFile >> word1 >> word2 >> word3;  //Position
		x = atof(word1.c_str());
		y = atof(word2.c_str());
		z = atof(word3.c_str());
		Point cam_p(x,y,z);

		sceneFile >> word1 >> word1 >> word2 >> word3; // Rotation
		x = (atof(word1.c_str()) * PI) * invRad;
		y = (atof(word2.c_str()) * PI) * invRad;
		z = (atof(word3.c_str()) * PI) * invRad;

		sceneFile >> word1 >> word1;
		aov = (atof(word1.c_str()) * PI) * invRad;

		sceneFile >> word1 >> word1 >> word2 >> word3; // Background color
		double r = atof(word1.c_str());
		double g = atof(word2.c_str());
		double b = atof(word3.c_str());

		bgcolor = Pixel(r,g,b);

		c.SetCam(aov,cam_p,x,y,z,r,g,b);

		// read in lights
		sceneFile >> word1 >> word1;
		numLights = atoi(word1.c_str());
		lights = new Light[numLights];

		for(int i = 0; i < numLights; i++)
		{
			sceneFile >> word1 >> word1; // "Light: Name"
			string name =  word1;
			sceneFile >> word1 >> word1; // Light Type
			int type = atoi(word1.c_str());

			sceneFile >> word1 >> word1 >> word2 >> word3; // Rotation
			x = (atof(word1.c_str()) * PI) * invRad;
			y = (atof(word2.c_str()) * PI) * invRad;
			z = (atof(word3.c_str()) * PI) * invRad;

			double tx = 0.;
			double ty = sin(x);
			double tz = -cos(x);

			double nx = tz*sin(y) + tx*cos(y);
			double ny = ty;
			double nz = tz*cos(y) - tx*sin(y);

			tx = nx*cos(z) - ny*sin(z);
			ty = nx*sin(z) + ny*cos(z);
			tz = nz;

			Vector d(tx, ty, tz);

			sceneFile >> word1 >> word1 >> word2 >> word3; // Color
			double r,g,b;
			r = atof(word1.c_str());
			g = atof(word2.c_str());
			b = atof(word3.c_str());

			sceneFile >> word1 >> word1;  // Light Intensity
			double inten = atof(word1.c_str());

			sceneFile >> word1 >> word1 >> word2 >> word2 >> word3 >> word3;  // Diff and Spec
			bool diff = bool(atoi(word1.c_str()));
			bool spec = bool(atoi(word2.c_str()));
			bool shad = bool(atoi(word3.c_str()));

			sceneFile >> word1 >> word1 >> word2 >> word2 >> word3 >> word3;
			double shadAngle = atof(word1.c_str());
			int shadRays = atoi(word2.c_str());
			int shadDep = atoi(word3.c_str());

			// Point light info
			tx = ty = tz = 0.;
			int dec = 0;
			double cone = 0.;
			double drop = 0.;
			double fullAngle = 0.;
			if(type > 0)
			{
				sceneFile >> word1 >> word1 >> word2 >> word3;  // Position
				tx = atof(word1.c_str());
				ty = atof(word2.c_str());
				tz = atof(word3.c_str());

				sceneFile >> word1 >> word1; // Decay Rate
				dec = atoi(word1.c_str());

				// Spot light info
				if(type > 1)
				{
					sceneFile >> word1 >> word1 >> word2 >> word2 >> word3 >> word3;
					//cone = cos(atof(word1.c_str()) * PI) / 180.;
					cone = atof(word1.c_str()) * 0.5f;
					double pen = atof(word2.c_str());
					drop = atof(word3.c_str());

					// ensure that the larger angle is stored in fullAngle
					if(pen > 0)
					{
						fullAngle = cos(((cone + pen) * PI) * invRad);
						cone = cos((cone * PI) * invRad);
					}
					else
					{
						fullAngle = cos((cone * PI) * invRad);
						cone = cos(((cone + pen) * PI) * invRad);
					}
				}
			}

			Point p(tx,ty,tz);

			Light lig(name, p, d, r, g, b, inten, type, diff, spec, shad, shadAngle, shadRays, shadDep, dec, cone, fullAngle, drop);
			lights[i] = lig;
		}

		// read in Textures
		sceneFile >> word1 >> word1;
		numTex = atoi(word1.c_str());
		tex = new Texture[numTex];

		for(int i = 0; i < numTex; i++)
		{
			sceneFile >> word1;
			tex[i].filename = word1;
			if(!tex[i].ReadFile())
				return false;
		}

		// read in Materials
		sceneFile >> word1 >> word1;
		numMats = atoi(word1.c_str());
		mats = new Material[numMats];

		sceneFile >> word1;
		for(int i = 0; i < numMats; i++)
		{
			sceneFile >> word1 >> word2 >> word2;  // Name and Type
			string name = word1;
			bool type = bool(atoi(word2.c_str()));

			sceneFile >> word1 >> word1; // Diffuse
			double diff = atof(word1.c_str());

			int diffTex = -1;
			int opacTex = -1;
			int ambiTex = -1;
			int emitTex = -1;
			int specTex = -1;
			int bumpTex = -1;
			double freqX =  0.;
			double freqY =  0.;
			double freqZ =  0.;
			double offX = 0.;
			double offY = 0.;
			double offZ = 0.;
			double bumpDepth = 0.;

			double dr, dg, db;
			dr = dg = db = 0.;
			sceneFile >> word1 >> word1; // Diff Color
			if(word1 == "Texture")
			{
				sceneFile >> word1;
				diffTex = atoi(word1.c_str());
			}
			else
			{
				sceneFile >> word2 >> word3;
				dr = atof(word1.c_str());
				dg = atof(word2.c_str());
				db = atof(word3.c_str());
			}

			double tr, tg, tb;
			tr = tg = tb = 0.;
			sceneFile >> word1 >> word1;  // Transparency
			if(word1 == "Texture")
			{
				sceneFile >> word1;
				opacTex = atoi(word1.c_str());
			}
			else
			{
				sceneFile >> word2 >> word3;
				tr = atof(word1.c_str());
				tg = atof(word2.c_str());
				tb = atof(word3.c_str());
			}

			double ar, ag, ab;
			ar = ag = ab = 0.;
			sceneFile >> word1 >> word1;  // Ambient Color
			if(word1 == "Texture")
			{
				sceneFile >> word1;
				ambiTex = atoi(word1.c_str());
			}
			else
			{
				sceneFile >> word2 >> word3;
				ar = atof(word1.c_str());
				ag = atof(word2.c_str());
				ab = atof(word3.c_str());
			}

			double er, eg, eb;
			er = eg = eb = 0.;
			sceneFile >> word1 >> word1;  // Incandescence
			if(word1 == "Texture")
			{
				sceneFile >> word1;
				emitTex = atoi(word1.c_str());
			}
			else
			{
				sceneFile >> word2 >> word3;
				er = atof(word1.c_str());
				eg = atof(word2.c_str());
				eb = atof(word3.c_str());
			}

			sceneFile >> word1 >> word1 >> word2 >> word2; //Refractions
			bool refr = bool(atoi(word1.c_str()));
			double refInd = 1. / atof(word2.c_str());

			// Specular Info
			double spec, sr, sg, sb, refl;
			spec = sr = sg = sb = refl = 0.;
			if(type != 0)
			{
				sceneFile >> word1 >> word1; // Phong Power
				//spec = atof(word1.c_str());
				//spec = sin(( max(min(atof(word1.c_str()),90.),0.) * PI) * invRad);
				spec = (( max(min(atof(word1.c_str()),90.),0.) * PI) * invRad);

				sceneFile >> word1 >> word1;  // Spec Color
				if(word1 == "Texture")
				{
					sceneFile >> word1;
					specTex = atoi(word1.c_str());
				}
				else
				{
					sceneFile >> word2 >> word3;
					sr = atof(word1.c_str());
					sg = atof(word2.c_str());
					sb = atof(word3.c_str());
				}

				sceneFile >> word1 >> word1;  // Reflectivity
				refl = atof(word1.c_str());
			}

			sceneFile >> word1;
			if(word1 == "Bump")
			{
				sceneFile >> word2 >> word3 >> word1;
				bumpTex = atoi(word2.c_str());
				bumpDepth = atof(word3.c_str());
			}
			else if(word1 == "Noise")
			{
				sceneFile >> word1 >> word2 >> word3;
				freqX = atof(word1.c_str());
				freqY = atof(word2.c_str());
				freqZ = atof(word3.c_str());

				sceneFile >> word1 >> word2 >> word3;
				offX = atof(word1.c_str());
				offY = atof(word2.c_str());
				offZ = atof(word3.c_str());

				sceneFile >> word2 >> word1;
				bumpDepth = atof(word2.c_str());
			}

			// Constructor:      diffuse         transparen   ambient     specular
			// int dtex = -1, int stex = -1, int otex = -1, int atex = -1, int etex = -1, int btex = -1)
			Material mat(name, dr, dg, db, diff, tr, tg, tb, ar, ag, ab, er, eg, eb, refr, refInd,
				sr, sg, sb, spec, refl, type, diffTex, specTex, opacTex, ambiTex, emitTex, bumpTex,
				freqX, freqY, freqZ, offX, offY, offZ, bumpDepth);
			mats[i] = mat;
		}

		// read in shapes
		sceneFile >> word1;

		numShapes = atoi(word1.c_str());
		shapes = new Shape*[numShapes];
		emitters = new int[numShapes];
		numEmit = 0;

		sceneFile >> word1;

		for(int i = 0; i < numShapes; i++)
		{
			//sceneFile >> word1;

			if(word1 == "Mesh:")
			{
				sceneFile >> word1;  // Name
				string name = word1;

				sceneFile >> word1 >> word1;  // Tris
				int tris = atoi(word1.c_str());
				int numInds = tris*3;

				sceneFile >> word1;  // Indices
				int *indices = new int[numInds];
				for(int j = 0; j < numInds; j++)
				{
					sceneFile >> word1;
					indices[j] = atoi(word1.c_str());
				}

				sceneFile >> word1 >> word1;  // Vertices
				int numVerts = atoi(word1.c_str());
				Point *verts = new Point[numVerts];
				sceneFile >> word1;
				for(int j = 0; j < numVerts; j++)
				{
					sceneFile >> word1 >> word2 >> word3;
					x = atof(word1.c_str());
					y = atof(word2.c_str());
					z = atof(word3.c_str());

					Point vp(x,y,z);
					verts[j] = vp;
				}

				sceneFile >> word1;  // Normals
				Normal *norms = new Normal[numInds];
				for(int j = 0; j < numInds; j++)
				{
					sceneFile >> word1 >> word2 >> word3;
					x = atof(word1.c_str());
					y = atof(word2.c_str());
					z = atof(word3.c_str());

					Normal n(x,y,z);
					norms[j] = n;
				}

				sceneFile >> word1;  // UVs
				int numUVs = numInds*2;
				double *uvs = new double[numUVs];
				for(int j = 0; j < numUVs; j++)
				{
					sceneFile >> word1;
					uvs[j] = atof(word1.c_str());
				}

				sceneFile >> word1 >> word2;  // double sided
				bool doubsided = bool(atoi(word2.c_str()));

				sceneFile >> word1;
				bool primary = true;
				bool shadows = true;
				if(word1 == "Visible")
				{
					sceneFile >> word2;  // primary visibility
					primary = bool(atoi(word2.c_str()));

					sceneFile >> word1 >> word2;  // casts shadows
					shadows = bool(atoi(word2.c_str()));

					sceneFile >> word1;
				}

				sceneFile >> word1;  // Material ID
				int matId = -1;
				int j = 0;
				while(j < numMats && matId == -1)
				{
					if(mats[j].name == word1)
						matId = j;

					j++;
				}

				Mesh tMesh(tris, indices, verts, norms, uvs, numInds, numVerts, doubsided);
				tMesh.matId = matId;
				tMesh.primary = primary;
				tMesh.shadows = shadows;
				shapes[i] = tMesh.clone();

				if(mats[matId].er > 0 || mats[matId].eg > 0 || mats[matId].eb > 0 || mats[matId].emitTex != -1)
					emitters[numEmit++] = i;

				sceneFile >> word1;
			}

			else if(word1 == "Sphere:")
			{
			}
		}

		sceneFile.close();

		return true;
	}

	return false;
}
///////////////// END LOAD SCENE ///////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////


#endif
