import numpy as np
import math
import random
class Ball:
    
    def __init__(self, mass, radius, position,velocity,size):
        self.mass = mass
        self.radius = radius
        self.position = np.array(position)
        self.fixed_velocity=np.array([0,0])
        self.velocity = np.array([0,0])
        self.acceleration=np.array([0,0])
        self.random_point=np.array([0,0])
        self.fdr=False  # check if its inside the circle 

    def compute_fixed_velocity(self,position,velocity,size):
        position=self.position
        if not(self.fdr):#checking the direction of velocity to be assigned 
            if(position[1]<size/2):
                theta=math.atan((size/2-position[1])/(size-position[0]))
            else:
                theta=-math.atan((position[1]-size/2)/(size-position[0]))
            return np.array([math.cos(theta)*velocity,math.sin(theta)*velocity])
        else:    # code for velocity calculation of return 
            random_position=self.random_point[1]
            theta=math.atan((position[1]-random_position)/position[0])
            return np.array([-math.cos(theta)*velocity,-math.sin(theta)*velocity])# negative because direction is reversed 


    
    def compute_wall_acceleration(self,step,toug,A,B,k1,k2,size):
        position=self.position
        self.acceleration=np.add(self.acceleration,self.compute_wall(toug,A,B,k1,k2,position[0]%size,np.array([1,0])))
        self.acceleration=np.add(self.acceleration,self.compute_wall(toug,A,B,k1,k2,(size-position[1])%size,np.array([0,1])))
        self.acceleration=np.add(self.acceleration,self.compute_wall(toug,A,B,k1,k2,(size-position[0])%size,np.array([-1,0])))
        self.acceleration=np.add(self.acceleration,self.compute_wall(toug,A,B,k1,k2,position[1]%size,np.array([0,-1])))
        

    def compute_wall(self,toug,A,B,k1,k2,d,nij):
        r,v,m=self.radius,self.velocity,self.mass
        position=self.position
        ans=np.array([0,0])
        if(r>=d):            
            tij=np.array([-nij[1],nij[0]])
            delta_v=np.dot(v,tij)
            ans=np.add((k1/m)*nij+(k2/m)*delta_v*tij,ans,casting="unsafe")
            ans=np.add(ans,A/m*np.exp((r-d)/B)*nij,casting="unsafe")
            
        else:
            ans=np.add(ans,A/m*np.exp((r-d)/B)*nij,casting="unsafe")
        return ans



    def compute_acceleration(self,ball2,toug,A,B,k1,k2,size):
        self.fixed_velocity=self.compute_fixed_velocity(self.position,0.5,size)
        self.acceleration+=(self.fixed_velocity-self.velocity)/toug
        x1,x2=self.position,ball2.position
        r1,r2=self.radius,ball2.radius
        v1,v2=self.velocity,ball2.velocity
        d=np.linalg.norm(x1-x2)
        nij=(x1-x2)/d
        if(r1+r2>d):
            tij=np.array([-nij[1],nij[0]])
            delta_v=np.dot((v2-v1),tij)
            self.acceleration+=k1/self.mass*nij+k2/self.mass*delta_v*tij
            self.acceleration=np.add(A/self.mass*np.exp((r1+r2-d)/B)*nij,self.acceleration,casting="unsafe")
        else:
            self.acceleration=np.add(A/self.mass*np.exp((r1+r2-d)/B)*nij,self.acceleration,casting="unsafe")


    def compute_step(self, step,size):
        """Compute position of next step."""
        self.position = np.add(self.position,((step * self.velocity)+(1/2*self.acceleration*step*step)),casting="unsafe")
        self.velocity=np.add(self.velocity,step*self.acceleration)
        
    def  check_inside_fdr(self,size,radius):
        position=self.position
        po=np.array([size,size/2])
        if np.linalg.norm(position-po) < radius and not self.fdr:
            self.fdr=True
            self.random_point[1]+=random.randint(0,499)# assigning a random point on the 1st wall 
            print(self.random_point)#printing the position they are assigned to go

def solve_step_2(ball_list,step,size,obstacle_size,obstacle_radius):
    A=2000
    B=0.08
    k1=120000
    k2=240000
    toug=0.5
    for ball1 in ball_list:
        if ball1.mass==obstacle_size and ball1.radius==obstacle_radius:
            continue
        ball1.acceleration=0
        ball1.check_inside_fdr(size,20) # checking if ball inside the radius
        for ball2 in ball_list:
            if ball1 is not ball2:
                ball1.compute_acceleration(ball2,toug,A,B,k1,k2,size)
        ball1.compute_wall_acceleration(step,toug,A,B,k1,k2,size)
    for ball1 in ball_list:
        if ball1.mass==obstacle_size and ball1.radius==obstacle_radius:
            continue
        ball1.compute_step(step,size)
