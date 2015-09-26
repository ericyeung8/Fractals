#include <math.h>
#include "assignment3.h"
#include "init.h"
#include <queue>
#include <GL/glut.h>

//Global Variables
vector<Pt> drawpoints;
vector<Matrix> Trans;
queue<Pt> PtSet;
int color = 0;

//Calculate the inverse of matrix 3x3
Matrix inverse(Matrix m)
{
	Matrix rvalue;
	float para;

	para = 1.0 / (m.data[0][0] * m.data[1][1] - m.data[0][0] * m.data[1][2] - m.data[0][1] * m.data[1][0] + m.data[0][1] * m.data[1][2] + m.data[0][2] * m.data[1][0] - m.data[0][2] * m.data[1][1]);
	rvalue.data[0][0] = para*(m.data[1][1] - m.data[1][2]);
	rvalue.data[0][1] = para*(m.data[0][2] - m.data[0][1]);
	rvalue.data[0][2] = para*(m.data[0][1] * m.data[1][2] - m.data[0][2] * m.data[1][1]);
	rvalue.data[1][0] = para*(m.data[1][2] - m.data[1][0]);
	rvalue.data[1][1] = para*(m.data[0][0] - m.data[0][2]);
	rvalue.data[1][2] = para*(m.data[0][2] * m.data[1][0] - m.data[0][0] * m.data[1][2]);
	rvalue.data[2][0] = para*(m.data[1][0] - m.data[1][1]);
	rvalue.data[2][1] = para*(m.data[0][1] - m.data[0][0]);
	rvalue.data[2][2] = para*(m.data[0][0] * m.data[1][1] - m.data[0][1] * m.data[1][0]);

	return rvalue;
}

Matrix translate ( Vec v )
{
	Matrix rvalue;
	int i;

	for (i = 0; i <= 2; i++)
		rvalue.data[i][i] = 1;
	rvalue.data[0][2] = v.x;
	rvalue.data[1][2] = v.y;

	return rvalue;
}

Matrix rotate ( Pt p, float theta )
{
	Matrix rvalue;
	int i;

	for (i = 0; i <= 1; i++)
		rvalue.data[i][i] = cos(theta);
	rvalue.data[0][1] = -sin(theta);
	rvalue.data[0][2] = p.x + p.y*sin(theta) - p.x*cos(theta);
	rvalue.data[1][0] = sin(theta);
	rvalue.data[1][2] = p.y - p.y*cos(theta) - p.x*sin(theta);
	rvalue.data[2][2] = 1;

	return rvalue;
}

Matrix scale ( Pt p, float alpha )
{
	Matrix rvalue;
	int i;

	for (i = 0; i <= 1; i++)
		rvalue.data[i][i] = alpha;
	rvalue.data[0][2] = (1 - alpha)*p.x;
	rvalue.data[1][2] = (1 - alpha)*p.y;
	rvalue.data[2][2] = 1;

	return rvalue;
}

Matrix nscale ( Pt p, Vec v, float alpha )
{
	Matrix rvalue;

	rvalue.data[0][0] = (alpha - 1)*v.x*v.x + 1;
	rvalue.data[0][1] = (alpha - 1)*v.x*v.y;
	rvalue.data[0][2] = v.x*(p.x*v.x + p.y*v.y)*(1 - alpha);
	rvalue.data[1][0] = (alpha - 1)*v.x*v.y;
	rvalue.data[1][1] = (alpha - 1)*v.y*v.y + 1;
	rvalue.data[1][2] = v.y*(p.x*v.x + p.y*v.y)*(1 - alpha);
	rvalue.data[2][2] = 1;

	return rvalue;
}

Matrix image ( Pt p1, Pt p2, Pt p3, Pt q1, Pt q2, Pt q3 )
{
	Matrix rvalue;
	Matrix q, p;
	int i;

	p.data[0][0] = p1.x;
	p.data[0][1] = p2.x;
	p.data[0][2] = p3.x;
	p.data[1][0] = p1.y;
	p.data[1][1] = p2.y;
	p.data[1][2] = p3.y;

	q.data[0][0] = q1.x;
	q.data[0][1] = q2.x;
	q.data[0][2] = q3.x;
	q.data[1][0] = q1.y;
	q.data[1][1] = q2.y;
	q.data[1][2] = q3.y;

	for (i = 0; i <= 2; i++)
	{
		p.data[2][i] = 1;
		q.data[2][i] = 1;
	}

	rvalue = compose(q, inverse(p));

	return rvalue;
}

Matrix compose ( Matrix m2, Matrix m1 )
{
	Matrix rvalue;
	int row, col;
	int i;

	for (row = 0; row <= 2; row++)
	for (col = 0; col <= 2; col++)
	for (i = 0; i <= 2; i++)
		rvalue.data[row][col] += m2.data[row][i] * m1.data[i][col];

	return rvalue;
}

//Calculate the product of matrix 3x3 and matrix 3x1
Points Multiply(Matrix m, Pt p)
{
	Points rvalue;

	rvalue.x = m.data[0][0] * p.x + m.data[0][1] * p.y + m.data[0][2];
	rvalue.y = m.data[1][0] * p.x + m.data[1][1] * p.y + m.data[1][2];
	rvalue.z = 1;

	return rvalue;
}

//Store the Condensation Set
void setCondensationSet ( vector<Pt> pts )
{
	unsigned int i;

	for (i = 0; i < pts.size(); i++)
	{
		drawpoints.push_back(pts.back());
		pts.insert(pts.begin(), pts.back());
		pts.pop_back();
	}

	while (!pts.empty())
	{
		PtSet.push(pts.back());
		pts.pop_back();
	}

}

//Do the recursion with on condensation set
void Iteration( vector<Matrix> transformations, Pt p, int n)
{
	unsigned int size, i;

	size = transformations.size();

	Points pts;

	if (n > 0)
	{
		for (i = 0; i <= size-1; i++)	//for each transformation
		{
			pts = Multiply(transformations[i], p);
			Pt temp(pts.x, pts.y);
			Iteration(transformations, temp, n-1);
		}
	}
	else
		drawpoints.push_back(p);
}

//Do the recursion with condensation set
void IterationSet(vector<Matrix> transformations, vector<Pt> p, int n)
{
	unsigned int size, setsize, i, j;
	vector<Pt> newp;
	vector<Pt> q;
	Points pts;

	size = transformations.size();
	setsize = PtSet.size();

	if (n > 0)
	{
		for (i = 0; i <= size - 1; i++)	//for each transformation
		{
			for (j = 0; j <= p.size() - 1; j++)	//for each point in the condensation set
			{
				q.push_back(p.back());
				p.insert(p.begin(), p.back());
				p.pop_back();
			}

			while (!q.empty())	//store each point in the tree and will display all
			{
				pts = Multiply(transformations[i], q.back());
				q.pop_back();
				Pt temp(pts.x, pts.y);
				drawpoints.push_back(temp);
				newp.push_back(temp);
			}
			IterationSet(transformations, newp, n - 1);
		}
		//for (i = 0; i <= size - 1; i++)
		//{
		//	for (j = 0; j <= setsize - 1; j++)
		//	{
		//		pts = Multiply(transformations[i], PtSet.front());
		//		PtSet.push(PtSet.front());
		//		PtSet.pop();
		//		Pt temp(pts.x, pts.y);
		//		newp.push_back(temp);
		//	}
		//	IterationSet(transformations, newp, n - 1);
		//}
	}
	//else
	//{
	//	while (!p.empty())
	//	{
	//		drawpoints.push_back(p.back());
	//		p.pop_back();
	//	}
	//}
}

//Store the transformations
void setIATTransformations ( vector<Matrix> transformations )
{
	while (!transformations.empty())
	{
		Trans.push_back(transformations.back());
		transformations.pop_back();
	}
}

//The entrance of recursion
void ExecuteIteration(void)
{
	int n;
	unsigned int i;
	Pt p;

	if (Trans.size() <= 3)
		n = 10;
	else
		n = 6;
	if (color == 5)
		n = 30;
	else if (color == 6)
		n = 38;

	if (PtSet.empty() && !Trans.empty())
		Iteration(Trans, p, n);
	else if (!Trans.empty())
	{
		vector<Pt> pts;
		for (i = 0; i <= PtSet.size() - 1; i++)
		{
			pts.push_back(PtSet.front());
			PtSet.push(PtSet.front());
			PtSet.pop();
		}
		IterationSet(Trans, pts, n-2);
	}

	color = 0;
}

//Six fractals
void Figure1()
{
	/*figure 1*/
	vector<Matrix> iat;
	iat.push_back(scale(Pt(-0.9, -0.9), 0.5));
	iat.push_back(scale(Pt(0.9, -0.9), 0.5));
	iat.push_back(compose(rotate(Pt(0, 0), 3.1415926 / 2), scale(Pt(0.9, -0.9), 0.5)));

	setIATTransformations(iat);

	vector<Pt> pts;
	// no condensation set
	setCondensationSet(pts);

}

void Figure2()
{
	/*figure 2*/
	vector<Matrix> iat;
	iat.push_back(scale(Pt(-0.9, 0.9), 0.5));
	iat.push_back(scale(Pt(0.9, -0.9), 0.5));

	setIATTransformations(iat);

	vector<Pt> pts;
	pts.push_back(Pt(-0.9, -0.9));
	pts.push_back(Pt(-0.9, 0));
	pts.push_back(Pt(0, 0));
	pts.push_back(Pt(0, -0.9));
	
	setCondensationSet(pts);

}

void Figure3()
{
	/*figure 3*/
	vector<Matrix> iat;
	float asin1, asin2, acos1, acos2, s;

	asin1 = 0.9 * sin(18.0 / 180 * 3.141592653);
	asin2 = 0.9 * sin(54.0 / 180 * 3.141592653);
	acos1 = 0.9 * cos(18.0 / 180 * 3.141592653);
	acos2 = 0.9 * cos(54.0 / 180 * 3.141592653);
	s = 1 / (cos(72.0 / 180 * 3.141592653)*2 + 2);

	iat.push_back(scale(Pt(-acos1,  asin1), s));
	iat.push_back(scale(Pt( 0.0,  0.9), s));
	iat.push_back(scale(Pt( acos1,  asin1), s));
	iat.push_back(scale(Pt( acos2, -asin2), s));
	iat.push_back(scale(Pt(-acos2, -asin2), s));
	iat.push_back(compose(scale(Pt(0,0),s),rotate(Pt(0,0),3.141592653)));

	setIATTransformations(iat);

	vector<Pt> pts;
	// no condensation set
	setCondensationSet(pts);

}

void Figure4()
{
	/*figure 4*/
	vector<Matrix> iat;
	iat.push_back(scale(Pt(-0.9,   0.0), 1.0 / 3.0));
	iat.push_back(scale(Pt(-0.45,  0.9), 1.0 / 3.0));
	iat.push_back(scale(Pt( 0.45,  0.9), 1.0 / 3.0));
	iat.push_back(scale(Pt( 0.9,   0.0), 1.0 / 3.0));
	iat.push_back(scale(Pt( 0.45, -0.9), 1.0 / 3.0));
	iat.push_back(scale(Pt(-0.45, -0.9), 1.0 / 3.0));

	setIATTransformations(iat);

	vector<Pt> pts;
	// no condensation set
	setCondensationSet(pts);

}

void Figure5()
{
	vector<Matrix> iat;
	iat.push_back(compose(rotate(Pt(0.0, 0.0), 11.25 / 180 * 3.1415926), scale(Pt(0.0, 0.0), 2.5 / 3)));

	setIATTransformations ( iat );

	vector<Pt> pts;
	//condensation set
	pts.push_back(Pt(-0.3, 0.8));
	pts.push_back(Pt(0, 0.9));
	pts.push_back(Pt(0.3, 0.8));
	pts.push_back(Pt(0.6*sin(45.0 / 180 * 3.1415926), 0.6*sin(45.0 / 180 * 3.1415926)));
	pts.push_back(Pt(0.8, 0.3));
	pts.push_back(Pt(0.9, 0));
	pts.push_back(Pt(0.8, -0.3));
	pts.push_back(Pt(0.6*sin(45.0 / 180 * 3.1415926), -0.6*sin(45.0 / 180 * 3.1415926)));
	pts.push_back(Pt(0.3, -0.8));
	pts.push_back(Pt(0, -0.9));
	pts.push_back(Pt(-0.3, -0.8));
	pts.push_back(Pt(-0.6*sin(45.0 / 180 * 3.1415926), -0.6*sin(45.0 / 180 * 3.1415926)));
	pts.push_back(Pt(-0.8, -0.3));
	pts.push_back(Pt(-0.9, 0));
	pts.push_back(Pt(-0.8, 0.3));
	pts.push_back(Pt(-0.6*sin(45.0 / 180 * 3.1415926), 0.6*sin(45.0 / 180 * 3.1415926)));
	setCondensationSet ( pts );

	color = 5;
}

void Figure6()
{
	vector<Matrix> iat;
	//iat.push_back(rotate(Pt(0.0, 0.0), 3.0/180*3.1415926));
	iat.push_back(image(Pt(-0.6, -0.6), Pt(0.6, -0.6), Pt(0.0, 0.6), Pt(-0.57, -0.63), Pt(0.63, -0.56), Pt(-0.03, 0.60)));

	setIATTransformations(iat);

	vector<Pt> pts;
	//condensation set
	pts.push_back(Pt(-0.6, -0.6));
	pts.push_back(Pt(0.6, -0.6));
	pts.push_back(Pt(0.0, 0.6));

	setCondensationSet(pts);

	color = 6;
}

// Draws the current IAT
void display ( void )
{
	Pt temp;
	unsigned int i;
	
	glColor3f(1, 0, 0);
	ExecuteIteration();

	if (!drawpoints.empty())
	{
		glClear(GL_COLOR_BUFFER_BIT);
		if (PtSet.empty() | PtSet.size() == 1)	//use glpoints
		{
			glBegin(GL_POINTS);
			while (!drawpoints.empty())
			{
				temp = drawpoints.back();
				drawpoints.pop_back();
				glVertex2f(temp.x, temp.y);

			}
			glEnd();
		}
		else
		{
			while (!drawpoints.empty())	//use gl_line_loop
			{
				glBegin(GL_LINE_LOOP);
				for (i = 0; i <= PtSet.size() - 1; i++)
				{
					temp = drawpoints.back();
					drawpoints.pop_back();
					glVertex2f(temp.x, temp.y);
				}
				glEnd();
			}
		}

		//clear all the global variables
		Trans.clear();
		while (!PtSet.empty())
			PtSet.pop();
	}

	glFlush ( );
}

/* do not modify the reshape function */
void reshape ( int width, int height )
{
	glViewport ( 0, 0, width, height );
	glMatrixMode ( GL_PROJECTION );
	glLoadIdentity ( );    
	gluOrtho2D (-1, 1, -1, 1);
	glMatrixMode ( GL_MODELVIEW );
    glLoadIdentity ( );
}

void keyboard(unsigned char key, int x, int y)
{
	switch (key)
	{
	case '1':
		Figure1();
		break;
	case '2':
		Figure2();
		break;
	case '3':
		Figure3();
		break;
	case '4':
		Figure4();
		break;
	case '5':
		Figure5();
		break;
	case '6':
		Figure6();
		break;
	}

	glutPostRedisplay();
}

int main ( int argc, char** argv )
{
	glutInit ( &argc, argv );
	glutInitDisplayMode ( GLUT_SINGLE | GLUT_RGB );
	glutInitWindowSize ( 500, 500 );
	glutInitWindowPosition ( 100, 100 );
	glutCreateWindow ( "Jiajun Yang - Homework 3" );
	init ( );
	glClear(GL_COLOR_BUFFER_BIT);
	glutDisplayFunc ( display );
	glutReshapeFunc ( reshape );
	glutKeyboardFunc( keyboard );
	glutMainLoop ( );
	return 0;
}
