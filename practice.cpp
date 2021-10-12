#include<iostream>
#include<algorithm>
#include<stack>
#include<math.h>
#include<cstdlib>
#include<set>
#include<stdlib.h>


using namespace std;

int n;
const double MAXNUM=1e10;

class c_point
{
public:
	double x,y;
	c_point(){};
	c_point(double x0, double y0){
		x = x0;
		y = y0;
	};
};

class polygon{
public:
	c_point *point_set;
	int num;
	int* line;
	polygon(c_point ps[], int* L, int n){
		point_set = ps;
		num = n;
		line = L;
	}
};

c_point p0;//最低点
c_point* points;//散点集
stack<c_point>convex_hull;//凸包

//显示栈
void check_stack(stack<c_point>points)
{
	cout<<"==========="<<endl;
	while(!points.empty()){
		printf("%2f %2f\n", points.top().x, points.top().y);
		points.pop();
	}
	cout<<"==========="<<endl;
}

//寻找p0
void find_p0()
{
	p0=points[0];
	int ii=0;
	for(int i=1;i<n;i++){
		if(points[i].y<p0.y){
			p0=points[i];
			ii=i;
		}else if(points[i].y==p0.y){
			if(points[i].x<p0.x){
				p0=points[i];
				ii=i;
			}
		}
	}
}

//极角排序
bool cmp(c_point p1,c_point p2)
{
	//p0排首位
	if(p1.x==p0.x&&p1.y==p0.y)return true;
	if(p2.x==p0.x&&p2.y==p0.y)return false;

	//计算极角（等于0则赋予一个极大值）
	double angle1=p1.x==p0.x?MAXNUM:(p1.y-p0.y)/(p1.x-p0.x);
	double angle2=p2.x==p0.x?MAXNUM:(p2.y-p0.y)/(p2.x-p0.x);
	//小于0则赋予一个更大的值
	if(angle1<0)angle1+=2*MAXNUM;
	if(angle2<0)angle2+=2*MAXNUM;
	
	//极角排序
	if(angle1<angle2)return true;
	else if(angle1==angle2){
		if(p1.y>p2.y)return true;
		else return false;
	}
	else return false;
}

//叉积
double cross(c_point p1, c_point p2, c_point p3)
{
	return (p2.x-p1.x)*(p3.y-p1.y)-(p3.x-p1.x)*(p2.y-p1.y);
}

//搜索凸包
void find_convex_hull()
{
	//p0和p1是凸包中的点
	convex_hull.push(points[0]);
	convex_hull.push(points[1]);

	int i=2;
	//p1,p2为栈顶两个节点
	c_point p1=points[0];
	c_point p2=points[1];
	while(i<n){
		//如果points[i]和points[i-1]在同一个角度，则不再对points[i]进行计算
		if((points[i-1].y-p0.y)*(points[i].x-p0.x)==(points[i-1].x-p0.x)*(points[i].y-p0.y)){
			i++;
			continue;
		}

		//如果叉积大于0，将当前点压入栈
		if (cross(p1, p2, points[i])>=0){
			//假设现在栈中为a,b,c,d,cross(c,d,e)大于等于0
			convex_hull.push(points[i]);//a,b,c,d,e,p1=c,p2=d
			p1=p2;//p1=d
			p2=convex_hull.top();//p2=e
			i++;
		}

		//如果叉积小于0，对栈中节点进行处理
		else{
			while(1){
				//假设现在栈中为a,b,c,d,cross(c,d,e)小于0
				convex_hull.pop();//a,b,c
				convex_hull.pop();//a,b
				p2=p1;//p2=c;
				p1=convex_hull.top();//p1=b
				convex_hull.push(p2);//a,b,c
				//cross(b,c,e)
				if(cross(p1,p2,points[i])>=0){
					convex_hull.push(points[i]);//a,b,c,e
					p1=p2;//p1=c
					p2=convex_hull.top();//p2=e
					i++;
					break;
				}
			}
		}
	}
}

// 判断点是否位于多边形内部
int is_inside(c_point p, polygon poly){
	int i,j,c=0;
	c_point *temp;
	for(i=0,j=poly.num-1; i<poly.num; j=i++){
		temp = poly.point_set;
		if(((temp[i].y>p.y)!=(temp[j].y>p.y))&&\
		(p.x<(temp[j].x-temp[i].x)*(p.y-temp[i].y)/(temp[j].y-temp[i].y)+temp[i].x))
		c = !c;
	}
	return c;
}

int main()
{
    // 输入数据，两组：两个多边形
	int n1;
    cin >> n1;
    c_point ps1[n1]; // 定义第一个多边形
	int L1[n1]; // 存储线段对应的端点
    for(int i=0;i < n1; i++){
        cin >> ps1[i].x >> ps1[i].y;
		L1[i] = n1-1-i;
    }
	polygon P1(ps1, L1, n1); // 创建多边形

    int n2;
    cin >> n2;
	c_point ps2[n2]; // 定义第二个多边形
	int L2[n2];
    for(int i=0; i < n2; i++){
        cin >> ps2[i].x >> ps2[i].y;
		L2[i] = n2-i-1;
    }
	polygon P2(ps2, L2, n2); // 创建多边形

    // 找到位于多边形内部的点
	stack<c_point> new_poly; // 收录重叠多边形的点
	for(int i=0; i<P1.num; i++){
		if(is_inside(P1.point_set[i], P2))
			new_poly.push(P1.point_set[i]);
	}

	for(int i=0; i<P2.num; i++){
		if(is_inside(P2.point_set[i], P1))
			new_poly.push(P2.point_set[i]);
	}

	// 找到两多边形的交点
	int i1, j1, i2, j2;
	for(i1=0, j1=P1.num-1;i1<P1.num;j1=i1++){
		for(i2=0, j2=P2.num-1; i2<P2.num; j2=i2++){
			// r1 r2是否异号
			double r1 = cross(P1.point_set[j1], P2.point_set[j2], P1.point_set[i1]);
			double r2 = cross(P1.point_set[j1], P2.point_set[i2], P1.point_set[i1]);
			// r3 r4是否异号
			double r3 = cross(P2.point_set[j2], P1.point_set[j1], P2.point_set[i2]);
			double r4 = cross(P2.point_set[j2], P1.point_set[i1], P2.point_set[i2]);

			double xj1, yj1, xi1, yi1;
			xj1 = P1.point_set[j1].x; yj1 = P1.point_set[j1].y; xi1 = P1.point_set[i1].x; yi1 = P1.point_set[i1].y;

			double xj2, yj2, xi2, yi2;
			xj2 = P2.point_set[j2].x; yj2 = P2.point_set[j2].y; xi2 = P2.point_set[i2].x; yi2 = P2.point_set[i2].y;

			if((r1 < 0 != r2 < 0) && (r3 < 0 != r4 < 0)){
				double a1, b1, c1, a2, b2, c2;
				a1 = P1.point_set[i1].y-P1.point_set[j1].y;
				b1 = P1.point_set[j1].x-P1.point_set[i1].x;
				c1 = P1.point_set[i1].x*P1.point_set[j1].y-P1.point_set[j1].x*P1.point_set[i1].y;

				a2 = P2.point_set[i2].y-P2.point_set[j2].y;
				b2 = P2.point_set[j2].x-P2.point_set[i2].x;
				c2 = P2.point_set[i2].x*P2.point_set[j2].y-P2.point_set[j2].x*P2.point_set[i2].y;

				double y = (c1*a2-c2*a1)/(a1*b2-a2*b1);
				double x = (c2*b1-c1*b2)/(a1*b2-a2*b1);

				new_poly.push(c_point(x, y));
			}
		}
	}

	n = new_poly.size();
	points = new c_point[n];
	int i=0;
	while (!new_poly.empty())
	{
		points[i++] = new_poly.top();
		new_poly.pop();
	}

    // 求凸包
    find_p0();//寻找p0
	sort(points,points+n,cmp);//按极角排序
	find_convex_hull();//搜索凸包

	// 根据凸包求面积
	double area = 0.;
	c_point p0 = convex_hull.top(); convex_hull.pop();
	c_point p1 = convex_hull.top(); convex_hull.pop();
	c_point p2;

	while(!convex_hull.empty()){
		p2 = convex_hull.top(); convex_hull.pop();
		area += 0.5*fabs(cross(p0, p1, p2));
		p1 = p2;
	}
	// check_stack(convex_hull);//显示结果
	printf("%f\n", area);
	cout << "finish" << endl;

    return 0;
}