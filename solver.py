import numpy as np
import math
class Ball:
    """Define physics of elastic collision."""
    
    def __init__(self, mass, radius, position,velocity,size):
        self.mass = mass
        self.radius = radius
        self.position = np.array(position)
        self.fixed_velocity=self.compute_fixed_velocity(position,velocity,size)
        self.velocity = np.array([0,0])
        self.acceleration=np.array([0,0])
        #self.vafter = np.copy(self.velocity) # temporary storage for velocity of next step

    def compute_fixed_velocity(self,position,velocity,size):
        position=self.position
        if(position[1]<size/2):
            theta=math.atan((size/2-position[1])/(size-position[0]))
        else:
            theta=-math.atan((position[1]-size/2)/(size-position[0]))
        return np.array([math.cos(theta)*velocity,math.sin(theta)*velocity])
    
    def compute_wall_acceleration(self,step,toug,A,B,k1,k2,size):
        position=self.position
        self.acceleration=np.add(self.acceleration,self.compute_wall(toug,A,B,k1,k2,position[0],np.array([1,0])))
        self.acceleration=np.add(self.acceleration,self.compute_wall(toug,A,B,k1,k2,size-position[1],np.array([0,1])))
        self.acceleration=np.add(self.acceleration,self.compute_wall(toug,A,B,k1,k2,size-position[0],np.array([-1,0])))
        self.acceleration=np.add(self.acceleration,self.compute_wall(toug,A,B,k1,k2,position[1],np.array([0,-1])))
        

    def compute_wall(self,toug,A,B,k1,k2,d,nij):
        r,v,m=self.radius,self.velocity,self.mass
        position=self.position
        ans=np.array([0,0])
        if(r>d):            
            tij=np.array([-nij[1],nij[0]])
            delta_v=np.dot(v,tij)
            ans=np.add((k1/self.mass)*nij+(k2/self.mass)*delta_v*tij,ans,casting="unsafe")
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
        d=np.linalg.norm(x2-x1)
        nij=(x1-x2)/d
        if(r1+r2>=d):
            tij=np.array([-nij[1],nij[0]])
            delta_v=np.dot((v2-v1),tij)
            self.acceleration+=k1/self.mass*nij+k2/self.mass*delta_v*tij
            self.acceleration=np.add(A/self.mass*np.exp((r1+r2-d)/B)*nij,self.acceleration,casting="unsafe")
        else:
            self.acceleration=np.add(A/self.mass*np.exp((r1+r2-d)/B)*nij,self.acceleration,casting="unsafe")


    def compute_step(self, step,size):
        """Compute position of next step."""
        #print(self.acceleration+"  "+self.position)
        self.position = np.add(self.position,step * self.velocity+(1/2*self.acceleration*(step)*step)%size,casting="unsafe")
        self.velocity=np.add(self.velocity,step*self.acceleration)
        


def solve_step_2(ball_list,step,size):
    A=2000
    B=0.08
    k1=120000
    k2=240000
    toug=0.5
    for ball1 in ball_list:
        ball1.acceleration=0
        ball1.compute_wall_acceleration(step,toug,A,B,k1,k2,size)
        for ball2 in ball_list:
            if ball1 is not ball2:
                ball1.compute_acceleration(ball2,toug,A,B,k1,k2,size)
    for ball1 in ball_list:
        ball1.compute_step(step,size)
