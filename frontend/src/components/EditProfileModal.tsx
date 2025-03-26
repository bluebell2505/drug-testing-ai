import React, { useState } from 'react';
import { 
  Dialog, 
  DialogContent, 
  DialogHeader, 
  DialogTitle,
  DialogFooter,
  DialogDescription
} from '@/components/ui/dialog';
import { Button } from '@/components/ui';
import { Input } from '@/components/ui/input';
import { Label } from '@/components/ui/label';
import { Textarea } from '@/components/ui/textarea';
import { useToast } from '@/hooks/use-toast';

interface ProfileData {
  user: {
    name: string;
    email: string;
    institution: string;
    role: string;
    joinedDate: string;
    initials: string;
  };
  preferences: {
    targetProteins: string[];
    optimizationGoals: string[];
    notificationSettings: string;
  };
  customization: {
    bio: string;
  };
  // Other profile data fields are omitted for simplicity
}

interface EditProfileModalProps {
  isOpen: boolean;
  onClose: () => void;
  profileData: ProfileData;
}

export const EditProfileModal: React.FC<EditProfileModalProps> = ({
  isOpen,
  onClose,
  profileData
}) => {
  const { toast } = useToast();
  const [formData, setFormData] = useState({
    name: profileData.user.name,
    institution: profileData.user.institution,
    role: profileData.user.role,
    email: profileData.user.email,
    bio: profileData.customization.bio,
    targetProteins: profileData.preferences.targetProteins.join(', '),
    optimizationGoals: profileData.preferences.optimizationGoals.join(', ')
  });

  const handleChange = (e: React.ChangeEvent<HTMLInputElement | HTMLTextAreaElement>) => {
    const { name, value } = e.target;
    setFormData(prev => ({ ...prev, [name]: value }));
  };

  const handleSubmit = (e: React.FormEvent) => {
    e.preventDefault();
    // In a real app, this would update the user profile in the database
    console.log('Updated profile data:', formData);
    
    // Show success message
    toast({
      title: "Success",
      description: "Profile updated successfully."
    });
    
    // Close the modal
    onClose();
  };

  return (
    <Dialog open={isOpen} onOpenChange={onClose}>
      <DialogContent className="sm:max-w-[425px] bg-gray-800 text-white border-gray-700">
        <DialogHeader>
          <DialogTitle className="text-xl font-semibold">Edit Profile</DialogTitle>
          <DialogDescription className="text-gray-400">
            Update your profile information.
          </DialogDescription>
        </DialogHeader>
        
        <form onSubmit={handleSubmit} className="space-y-4 pt-4">
          <div className="space-y-2">
            <Label htmlFor="name" className="text-white">Full Name</Label>
            <Input
              id="name"
              name="name"
              value={formData.name}
              onChange={handleChange}
              className="bg-gray-700 border-gray-600 text-white"
            />
          </div>
          
          <div className="space-y-2">
            <Label htmlFor="email" className="text-white">Email</Label>
            <Input
              id="email"
              name="email"
              type="email"
              value={formData.email}
              onChange={handleChange}
              className="bg-gray-700 border-gray-600 text-white"
            />
          </div>
          
          <div className="space-y-2">
            <Label htmlFor="institution" className="text-white">Institution</Label>
            <Input
              id="institution"
              name="institution"
              value={formData.institution}
              onChange={handleChange}
              className="bg-gray-700 border-gray-600 text-white"
            />
          </div>
          
          <div className="space-y-2">
            <Label htmlFor="role" className="text-white">Role</Label>
            <Input
              id="role"
              name="role"
              value={formData.role}
              onChange={handleChange}
              className="bg-gray-700 border-gray-600 text-white"
            />
          </div>
          
          <div className="space-y-2">
            <Label htmlFor="bio" className="text-white">Bio</Label>
            <Textarea
              id="bio"
              name="bio"
              value={formData.bio}
              onChange={handleChange}
              className="bg-gray-700 border-gray-600 text-white resize-none h-20"
            />
          </div>
          
          <div className="space-y-2">
            <Label htmlFor="targetProteins" className="text-white">Preferred Target Proteins (comma separated)</Label>
            <Input
              id="targetProteins"
              name="targetProteins"
              value={formData.targetProteins}
              onChange={handleChange}
              className="bg-gray-700 border-gray-600 text-white"
            />
          </div>
          
          <div className="space-y-2">
            <Label htmlFor="optimizationGoals" className="text-white">Default Optimization Goals (comma separated)</Label>
            <Input
              id="optimizationGoals"
              name="optimizationGoals"
              value={formData.optimizationGoals}
              onChange={handleChange}
              className="bg-gray-700 border-gray-600 text-white"
            />
          </div>
          
          <DialogFooter className="pt-4">
            <Button variant="outline" className="border-gray-700 text-white" onClick={onClose}>
              Cancel
            </Button>
            <Button type="submit" className="bg-blue-500 hover:bg-blue-600 text-white">
              Save Changes
            </Button>
          </DialogFooter>
        </form>
      </DialogContent>
    </Dialog>
  );
};
